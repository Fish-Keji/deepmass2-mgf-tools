# DeepMASS2 MGF Tools

将 MZmine 导出的 MGF 转为 DeepMASS2 输入，通过原始 Feature ID 将注释结果对应回 MZmine。

## 添加哪些信息

本工具确保每条谱图显式包含 **COMPOUND_NAME 和 PRECURSOR_MZ**，并保留正确的 IONMODE。

| 字段 | 处理规则 |
| --- | --- |
| `COMPOUND_NAME` | 直接使用原始 `FEATURE_ID`，不加前缀、不重新编号；已有名称不同则替换 |
| `PRECURSOR_MZ` | 从 `PEPMASS` 的第一项复制前体 m/z；已有时检查是否一致 |
| `IONMODE` | 统一为 `positive` 或 `negative`；缺失时须明确指定 |
| 原有字段与峰数据 | 保留 `FEATURE_ID`、`FEATURE_FULL_ID`、`PEPMASS` 等，峰数据和谱图顺序不变 |

例如 `FEATURE_ID=16` 对应 `COMPOUND_NAME=16`，DeepMASS2 导出的 `16.csv` 可对应回 MZmine 的 Feature ID 16。

## 处理前后

以下为同一谱图的字段节选。

**MZmine 原始字段：**

```text
FEATURE_ID=16
PEPMASS=126.99724
IONMODE=POSITIVE
```

**处理后字段：**

```text
FEATURE_ID=16
PEPMASS=126.99724
IONMODE=positive
PRECURSOR_MZ=126.99724
COMPOUND_NAME=16
```

## 系统兼容性

转换脚本适用于 **Windows、macOS、Linux**，需 **Python 3.8+**，仅使用标准库，无需安装其他依赖或 DeepMASS2。

已在 Windows 实测；macOS、Linux 尚未实测。此兼容性说明针对转换脚本，DeepMASS2 本体需单独安装。

## 用法

解压后，在脚本目录运行：

```console
python prepare_mgf.py input.mgf output_deepmass2.mgf
```

macOS / Linux 若使用 `python3` 命令，将上述 `python` 改为 `python3`。

仅当原文件缺失极性时，按实际采集模式补充（负离子改为 `negative`）：

```console
python prepare_mgf.py input.mgf output_deepmass2.mgf --ionmode positive
```

输出 MGF 可直接加载到 DeepMASS2。缺失或重复 `FEATURE_ID`、质量字段冲突时会报错；已有输出文件不覆盖，请换一个文件名。


工具版本 1.1.0 · MIT License
