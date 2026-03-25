# CPLkit

`CPLkit.py` は Gaussian の TD-DFT 出力ファイルから transition electric dipole moment density cube, transition magnetic dipole moment density cube, および CPL 関連の集計データを生成する Python スクリプトです。このスクリプトは励起状態情報を読み取り、cube グリッド上で TEDM と TMDM の密度を再構成し、CPL 関連量を CSV ファイルとして出力します。数値計算には NumPy を用います。分子軌道 cube は Gaussian `cubegen` によって生成することもできますし、事前に生成した MO cube を直接利用することもできます。

## 主な機能

- Gaussian TD-DFT 出力から TEDM density cube を生成します。
- Gaussian TD-DFT 出力から TMDM density cube を生成します。
- 解析した励起状態に対する CPL 関連量を CSV として自動出力します。
- 各遷移について発光速度定数を `Radiative Rate Constant (ns^-1)` として CSV に自動出力します。
- `cubegen` を用いた checkpoint ベースの MO cube 生成と、事前生成済み `mo<MO>.cube` の再利用の両方に対応します。
- GaussView との互換性のために Gaussian cube のメタデータを保持します。
- 投影後の total cube に加えて Cartesian 成分 cube を任意で出力できます。

## 必要環境

### Python

- Python 3.10 以降を推奨します。

### Python パッケージ

以下で依存関係を導入できます。

```bash
pip install -r requirements.txt
```

### 外部ソフトウェア

次のいずれかのワークフローが必要です。

1. Gaussian `cubegen` が利用可能であり、`--chk` で checkpoint または formatted checkpoint を指定する方法
2. `--mocubes_dir` で事前生成済みの MO cube ディレクトリを指定する方法

## インストール

```bash
git clone https://github.com/yuuyanagata/CPLkit
cd CPLkit
pip install -r requirements.txt
```

## 入力

`--log` には Gaussian TD-DFT の出力ファイルを指定します。

cube を生成する場合には、追加で次のいずれかが必要です。

- `--chk`  
  `cubegen` を用いて MO cube を生成するための checkpoint または formatted checkpoint

- `--mocubes_dir`  
  `mo<MO>.cube` を含むディレクトリ

## 基本的な使い方

### 1. CPL の CSV のみを出力する

```bash
python CPLkit.py --log sample.out --cpl_only
```

### 2. `cubegen` を用いて励起状態 1 の TEDM と TMDM cube を生成する

```bash
python CPLkit.py --log sample.out --state 1 --chk sample.fchk
```

### 3. 事前生成済み MO cube から cube を生成する

```bash
python CPLkit.py --log sample.out --state 1 --mocubes_dir ./sample_S1_mocubes
```

### 4. Cartesian 成分 cube も同時に出力する

```bash
python CPLkit.py --log sample.out --state 1 --chk sample.fchk --keep_components
```

## CSV に出力される発光速度定数

CPL CSV には，各励起状態遷移について `Radiative Rate Constant (ns^-1)` 列が含まれます。この値は，Gaussian TD-DFT 出力に含まれる振動子強度と遷移波長から，次式によって見積もられます。

```math
k_{\mathrm r}(\mathrm{ns}^{-1}) = 6.67 \times 10^{4} \frac{f}{\lambda^{2}(\mathrm{nm})}
```

ここで，`f` は Gaussian TD-DFT が与える振動子強度であり，`λ` は同じ CSV 行に出力される nm 単位の遷移波長です。この式は，Einstein の自然放出速度を波長ベースの実用式として表したものです。

## 主なコマンドライン引数

- `--log`  
  Gaussian TD-DFT の出力ファイル

- `--state`  
  cube を生成する励起状態番号。`--cpl_only` を使う場合には不要です。

- `--chk`  
  `cubegen` 用の checkpoint または formatted checkpoint

- `--mocubes_dir`  
  `mo<MO>.cube` を含むディレクトリ

- `--cubegen`  
  `cubegen` 実行ファイルのパス。既定値は `cubegen` です。

- `--cubegen_npts`  
  `cubegen` に渡す第 1 引数

- `--cubegen_grid`  
  `cubegen` に渡す残りのグリッド引数。既定値は `-3 h` です。

- `--outdir`  
  出力ディレクトリ

- `--outprefix`  
  出力プレフィックス。既定値は `S<state>` です。

- `--overwrite_mo_cubes`  
  既存の MO cube キャッシュが存在しても再生成します。

- `--keep_components`  
  投影後の total cube に加えて `x`, `y`, `z` 成分 cube を出力します。

- `--dtype`  
  内部演算の浮動小数点精度。`float32` または `float64` を指定できます。

- `--coords`  
  座標処理モード。`auto`, `aligned`, `general` を指定できます。

- `--cpl_only`  
  CPL の CSV のみを出力し、cube 生成を省略します。

- `--no_cpl_csv`  
  CPL CSV の出力を無効化します。

- `--cpl_csv_path`  
  CPL CSV の明示的な出力パス

## 出力ファイル

選択したオプションに応じて、以下のようなファイルが生成されます。

- `<logstem>-CPL.csv`。各遷移の CPL 関連量に加えて `Radiative Rate Constant (ns^-1)` 列を含みます。
- `<logstem>_S1_TEDM_total.cube`
- `<logstem>_S1_TMDM_total.cube`
- `<logstem>_S1_TEDM_x.cube`
- `<logstem>_S1_TEDM_y.cube`
- `<logstem>_S1_TEDM_z.cube`
- `<logstem>_S1_TMDM_x.cube`
- `<logstem>_S1_TMDM_y.cube`
- `<logstem>_S1_TMDM_z.cube`

`cubegen` を用いる場合には、次のような MO cube のキャッシュディレクトリも作成されます。

- `<logstem>_S1_mocubes/`

## リポジトリ構成

```text
CPLkit.py
LICENSE
README.md
README_ja.md
requirements.txt
```

## 引用とクレジット

このスクリプトを再利用または改変する場合には、学術上の慣行として元のリポジトリへのクレジットを明記し、変更の有無を示すことが望まれます。MIT ライセンス自体は学術的引用を義務づけませんが、ソフトウェアの複製物またはその重要部分には著作権表示およびライセンス文を残す必要があります。

## ライセンス

本リポジトリは MIT License の下で公開されます。

詳細は [LICENSE](LICENSE) を参照してください。
