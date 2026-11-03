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
