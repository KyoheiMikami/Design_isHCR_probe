# Design isHCR probe

[Mikami _et al._, _Fly_, 2025](https://doi.org/10.1080/19336934.2024.2428499)を参照。\
ネッパジーン社の[ISH palette](https://nepagene.jp/en/products/fluorescent-stain-insituhcr/ishpalette)に対応している。

## 要件

- python3
- Biopython

## 使い方

1. `Pick_ProbeRegion.py`でprobeの結合サイトを検索する
2. 検索結果をBlastにかけて、オフターゲットのない結合サイトを選ぶ
3. `MakeProbe.py`でprobeを作成する

### `Pick_ProbeRegion.py`について

option：

|||
|--|--|
|-i, —input_file|入力ファイルを指定する、必須|
|-m, —gc_min|結合サイトのGCの最小値を指定する、デフォルトは`45`|
|-M, —gc_max|結合サイトのGCの最大値を指定する、デフォルトは`55`|
|-t, —target_name|プローブの名称を指定する、必須|
|-o, —output_directory|出力先を指定する、デフォルトはカレントディレクトリ|

output：

|||
|--|--|
|`Target_ProbeRegion.csv`|配列やGCを確認できる`.csv`|
|`Target_ProbeRegion.fasta`|Blast検索用のマルチファスタ|

### `MakeProbe.py`について

option：

|||
|--|--|
|-i, —input_file|入力ファイルを指定する、必須|
|-t, —target_name|プローブの名称を指定する、必須|
|-H, —hairpin|Short hairpin amplifierの蛍光色素を指定する、必須|
|-f, —output_format|出力ファイルの形式を指定する、デフォルトは`.csv`|
|-o, —output_directory|出力先を指定する、デフォルトはカレントディレクトリ|
