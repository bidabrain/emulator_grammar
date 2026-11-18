# z=0 GSMF 标定工作流 (dual-channel SN feedback + spin-jet AGN)

用 SWIFT-Emulator 标定 SWIFT-LRD 的 6 个 feedback 参数(Table 2,ColibrePressure
模型),目标是 25 Mpc 盒在 z=0 的 galaxy stellar mass function (GSMF)。

## 总览

```
 design_gsmf.py          SWIFT              pgalf            measure_gsmf.py        train + calibrate
┌──────────────┐   ┌──────────────┐   ┌──────────────┐   ┌──────────────────┐   ┌──────────────────┐
│ Latin 超立方 │ → │ 48 次模拟     │ → │ galaxy       │ → │ catalogue → GSMF │ → │ GP emulator →    │
│ 48×run_X.yml │   │ (每个 yml)   │   │ catalogue    │   │ emulator YAML    │   │ 对观测标定        │
└──────────────┘   └──────────────┘   └──────────────┘   └──────────────────┘   └──────────────────┘
   design_parameters.json ─────────────────────────────────────────────────────┘ (参数,log 空间)
```

标定的 6 个参数(其余固定):

| Table 2 符号 | yaml key                                             | 采样范围       | 采样方式 |
|-------------|------------------------------------------------------|---------------|---------|
| f_kin       | `EAGLEFeedback:SNII_f_kinetic`                       | [0, 1]        | 线性     |
| ΔT_SNe      | `EAGLEFeedback:SNII_delta_T_K`                       | 10^[6.5, 8.0] | **log** |
| f_E,min     | `EAGLEFeedback:SNII_energy_fraction_min`             | [0.05, 0.4]   | 线性     |
| f_E,max     | `EAGLEFeedback:SNII_energy_fraction_max`             | [0.4, 8.0]    | 线性     |
| P_pivot     | `EAGLEFeedback:SNII_energy_fraction_P_pivot_K_p_cm3` | 10^[2.5, 4.5] | **log** |
| M_seed      | `SPINJETAGN:subgrid_seed_mass_Msun`                  | 10^[4, 6]     | **log** |

> `f_E,min` / `f_E,max` 的区间故意不重叠(···0.4] 和 [0.4···),保证每个样本
> `f_E,min ≤ f_E,max`(ColibrePressure logistic 的要求)。

---

## 环境准备 (一次)

```bash
pip install swiftemulator          # 或在本 repo: pip install -e .
pip install "pyDOE==0.3.8"         # 注意:直接 pip install pyDOE 会装成错误的 'pydoe' 1.5.0
```

本 repo 目录本身就是 `swiftemulator` 源码。若没正式安装,运行时在命令前加
`PYTHONPATH=.`(下文命令里的 `PYTHONPATH=.` 均可在正式安装后省略)。

训练阶段还需要 `george`、`SALib`、`corner`、`emcee`。

---

## Stage 1 — 生成参数文件 (`design_gsmf.py`)

```bash
cd /Volumes/Elements/scratch/emulator
PYTHONPATH=. python3 design_gsmf.py [/path/to/base.yml]
```

- base 模板默认用 `./templet.yml`,也可用参数覆盖。模板里 `SNII_energy_fraction_model`
  必须是 `ColibrePressure`。
- 产出到 `design_gsmf/`:
  - `run_0.yml … run_47.yml` — 48 个 SWIFT 参数文件(base + 覆盖那 6 个 key)
  - `design_parameters.json` — 48 组参数值(**采样空间**,log 参数存的是 log10 值)
  - `design_gsmf_corner.png` — 采样分布图
- 要改采样数 / 范围,编辑 `design_gsmf.py` 顶部的 `number_of_simulations` 与
  `parameter_limits`。

> ⚠ 重复运行会**重新随机采样并覆盖** `run_*.yml` 和 `design_parameters.json`。
> 一旦开始用某批文件跑模拟,就别再跑这个脚本,否则参数和模拟对不上。

---

## Stage 2 — 跑 SWIFT (你在集群上做)

对 48 个 `run_X.yml` 各跑一次 SWIFT 到 z=0。每次模拟放在各自的输出目录,
z=0 的 snapshot 号对所有 run 相同(由 `output_list.txt` 决定)。

---

## Stage 3 — 跑 galaxy finder (pgalf)

对每个 run 的 z=0 snapshot 跑 pgalf,得到 galaxy catalogue。参考
`pgalf_swift_grammar/run_pgalf.sh` 和其中的 `*.sbatch`:

```bash
./run_pgalf.sh -i <sim_dir> -o <out_dir> -s <snap> -n 64 -N 4 -m "srun -n"
```

产出 `<out_dir>/FoF_Data/FoF.<NNNNN>/GALCATALOG.LIST.<NNNNN>`(二进制,
每个 galaxy 含 `mstar`,单位 Msun/h)。

> pgalf 的宇宙学是**编译期硬编码**(`#define REQUIRED_OMEGA_M/…/_H0`),换宇宙学要改
> C 源码重编译。本标定里 48 个 run 宇宙学相同,无需改动。

> **必须用 pgalf_swift_grammar ≥ `2570de6`(并重新编译 NewDD)。** SWIFT snapshot 的
> 长度/质量不带 h,h 只存在 `/Cosmology/h`;旧版 NewDD 读不到它,默认 h=1,导致
> catalogue 实际是 cMpc / Msun(FoF `box_size` = BoxSize,而不是 BoxSize×h)。
> 由此测出的 GSMF:φ 偏低约 0.48 dex、M* 偏高约 0.16 dex、孔径只有标称的约 0.69 倍。
> 用旧 NewDD 跑过的 run 需从 NewDD 起重跑(输出到新目录,不要在旧目录上用 `-r`),
> 再重新测 GSMF。详见 `pgalf_swift_grammar/NewDD/NOTE_SWIFT_h_units.md`。
> 注意 FoF 头里的 `hubble`(69.33)是编译期写死的,不能用来判断。

---

## Stage 4 — 测量 GSMF (`measure_gsmf.py`)

把每个 run 的 pgalf catalogue 转成 emulator 要的 GSMF YAML。

```bash
cd /Volumes/Elements/scratch/emulator
python3 measure_gsmf.py \
  --pgalf-repo    /Volumes/Elements/scratch/pgalf_swift_grammar \
  --root-template "/你的/pgalf输出/run_{id}" \
  --snapshot      <z=0 snapshot 号> \
  --design        design_gsmf/design_parameters.json \
  --out-dir       gsmf_measured
```

- `--root-template` 里的 `{id}` 会被替换成 run id("0"…"47"),按你的实际目录布局填。
- run id 从 `design_parameters.json` 读,**和 emulator 的参数 key 一一对应**。
- 还没被 pgalf 处理的 run 会列为 MISSING 并跳过,可增量重跑。
- `--aperture-kpc`(默认 `50`):M* 取 3D 50 proper-kpc 孔径内的恒星质量(COLIBRE/观测约定,
  从 `GALFIND.DATA` 读);`0` = catalogue 的 total bound M*。可一次给多个值
  (如 `--aperture-kpc 0 30 50 100`),每个 run 的 `GALFIND.DATA` 只读一遍,结果分别写到
  `gsmf_measured/ap_<X>kpc/`(0 → `ap_total/`);训练时选其中一个目录。
- `--pgalf-repo` 需 ≥ `3d0b545`(自动识别 `GALFIND.DATA` 记录长度;旧版写死长度,
  NCHEM/NDUST 编译下会报 `Truncated subhalo header`),多孔径需 ≥ `ba53e2d`;
  版本不够时脚本会直接报错提示 `git pull`。
- **单位自检**:每个 run 测量前,脚本把 FoF 的 `box_size` 和 `run_pgalf.sh` 链接的
  `<root>/snap_XXXX.hdf5` 里的 `BoxSize × h` 比对;若是旧 NewDD 产物(box_size = BoxSize)
  或两者都对不上,直接报错停止,不写 YAML。找不到该链接(或没有 `h5py`)的 run 会列为
  "little-h check SKIPPED";`--skip-h-check` 可关闭此检查。
- `--center com|galcenter`(默认 `com`):孔径中心。`com` = subhalo 所有粒子的质心;
  `galcenter` = GalCenter 输出的恒星中心(`GALFIND.CENTER`,需跑过 galcenter 阶段),
  缺失的 subhalo 回退到质心。`galcenter` 更接近观测做法,能避免质心偏离恒星主体时
  孔径漏掉恒星。
- **Zoom-in 模拟**:只有高分辨区域里有 galaxy,φ 必须按高分辨区域的体积归一化。
  用 `--volume-frac F`(φ 用 `F × 盒子体积`,如 1/4 × 1/4 × 1 的高分辨长条为 1/16;
  若扣掉边缘污染层则更小),否则 φ 会低 1/F 倍。`--max-lowres-frac F` 去掉低分辨 DM
  质量占比超过 F 的 galaxy(两种 finder 都支持;`0` = 只保留完全未污染的)。两者默认
  关闭,即均匀周期盒。
- **改用 HBT-HERONS**(`--finder hbt --hbt-repo <HBT-HERONS-GRAMMAR>`):Stage 3 改跑
  `HBT-HERONS-GRAMMAR/grammar/run_hbt.sh`(要求 SWIFT 每个 snapshot 都带 FoF 的
  `FOFGroupIDs`,且需从早期 snapshot 一路跑到 z=0,详见该 repo 的 README),
  `--root-template` 指向各 run 的 HBT 输出目录。YAML 格式和单位不变;galaxy = HBT 自束缚
  subhalo(≥20 个束缚粒子),M* 只算束缚恒星,孔径中心默认为最束缚粒子
  (`--center mostbound`,也可 `com`)。HBT 从 SWIFT 元数据读 h 和单位,不需要 little-h 自检。
  两种 finder 的结果请写到不同的 `--out-dir`。
- 产出 `gsmf_measured/run_X.yml`,结构:

```yaml
stellar_mass_function:
  lines:
    mass_function:
      centers:       [...]          # M* bin 中心 (线性 Msun)
      centers_units: Solar_Mass
      values:        [...]          # phi (线性 Mpc^-3 dex^-1)
      values_units:  "1/Mpc**3"
      scatter:       [...]          # phi 的 Poisson 误差 (线性)
    histogram:
      values:        [...]          # 每个 bin 的 galaxy 原始计数 (仅参考)
```

**单位**:pgalf 全程 h-单位;本脚本统一转成**物理单位**(除掉 h,Msun / Mpc⁻³),
以便下一步直接和观测(Baldry+2012 / Driver+2022,物理 Msun)对齐,不用再绕 h。

---

## Stage 5 — 训练 + 标定 (待写)

读入参数与 GSMF,训练 GP,再对观测 z=0 GSMF 做 penalty / MCMC:

```python
from swiftemulator.io.swift import load_parameter_files, load_pipeline_outputs
from swiftemulator.emulators.gaussian_process import GaussianProcessEmulator

# 参数:直接读采样空间的 design_parameters.json (或用 load_parameter_files 读 run_*.yml)
# GSMF:
values, units = load_pipeline_outputs(
    filenames={rid: f"gsmf_measured/run_{rid}.yml" for rid in run_ids},
    scaling_relations=["stellar_mass_function"],
    log_independent=["stellar_mass_function"],
    log_dependent=["stellar_mass_function"],
    sanitize_output=True,            # 丢掉 N 过小导致的非有限 bin
)
# emulator = GaussianProcessEmulator(); emulator.fit_model(spec, params, values[...])
# 然后用 swiftemulator.comparison.penalty 对观测 GSMF 做标定
```

> ⚠ **swiftemulator 的坑**:`load_pipeline_outputs` 若同时传 `histogram_for_error`
> 和 `log_dependent`,values 不会被取 log(`swift.py:155` 被跳过),得到"线性值 + log
> 误差"的不一致数据。**用文件里的 `scatter` 字段 + `sanitize_output=True`,不要用
> `histogram_for_error`。**

---

## 关键文件

| 文件 | 作用 |
|------|------|
| `design_gsmf.py`             | Stage 1:Latin 超立方 → 48 个 SWIFT 参数文件 |
| `measure_gsmf.py`            | Stage 4:pgalf catalogue → emulator GSMF YAML |
| `design_gsmf/`               | Stage 1 产物(run_*.yml, design_parameters.json, corner) |
| `gsmf_measured/`             | Stage 4 产物(每个 run 一个 GSMF YAML) |
| `templet.yml`                | SWIFT base 模板(ColibrePressure) |
