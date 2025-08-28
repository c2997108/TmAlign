# TmAlign

TmAlign — Local Duplex Stability Search

概要
- クエリ配列内の各 k-mer ウィンドウを使って FASTA の各エントリを走査し、近傍の二本鎖安定性（Tm）を評価します。
- 左右フランクで 1 塩基のバルジ（ギャップ）を許容し、固定ペナルティ（10°C）を適用します。
- FASTA のエントリ単位で並列化し、各エントリの計算が終わり次第すぐに結果を出力します。

必要要件
- Python 3.8+
- リポジトリ内の `Santalucia_NN_Tm.py`（サンタルチア最近接モデル）[ThermoAlign](https://github.com/drmaize/ThermoAlign)より

基本的な使い方
- 最小実行例:
  - `python3 TmAlign.py --query CGCCTGTTTATCAAAAACAT --fasta test.fa --k 4`
- 並列実行（CPU数に合わせ自動）:
  - `python3 TmAlign.py --query CGCCTGTTTATCAAAAACAT --fasta test.fa --k 4 --cpus 0`
- しきい値を調整（低い Tm も出力、同一率 60%以上）:
  - `python3 TmAlign.py --query CGCCTGTTTATCAAAAACAT --fasta test.fa --k 4 --min-tm 0 --min-identity 60 --cpus 4`

オプション
- `--fasta`: 入力 FASTA ファイルパス
- `--query`: クエリ配列
- `--k`: k-mer 長（デフォルト: 4）
- `--min-tm`: 出力する最小 Tm（°C、デフォルト: 30.0）
- `--min-identity`: 出力する最小一致率（% 、デフォルト: 50.0）
- `--cpus`: 並列プロセス数。`0`=CPU数、`1`=並列なし、`2+`=明示した数（デフォルト: 0）

出力形式（タブ区切り）
- 列: `Tm  contig  start  end  strand  identity(%)  query_align  db_align`
- `Tm`: 小数点以下 2 桁
- `start`/`end`: 1-based inclusive 座標（DB 側）
- `strand`: `+` はクエリ、`-` は `revcomp(query)` でのヒット
- `identity(%)`: 整列長を分母に、同一文字の位置数を分子として算出（ギャップ位置も分母に含む）
- `query_align` / `db_align`: バルジは相手側に `-` として表示されます

並列実行と並び順
- 並列単位は FASTA のエントリ（レコード）ごとです。
- 各エントリ内の結果は Tm の降順でソートして出力します。
- エントリ間では「完了順」に出力されます（入力順は保証しません）。

アルゴリズムの要点
- クエリの各 k-mer ウィンドウをキーに一致位置を探索（オーバーラップ含む）。
- 左右フランクについて、ベースライン（ギャップなし）・DB 側 1 塩基バルジ・クエリ側 1 塩基バルジを評価し、最も高い Tm を採用（バルジには 10°C のペナルティ）。
- 採用整列に対して `Santalucia_NN_Tm.py` のモデルで Tm を計算します。

ヒント・注意
- 出力はエントリの計算が完了した単位で順次流れます。出力を入力順に強制したい場合はご相談ください。
- 同一率の算出はギャップ位置も分母にカウントします（`a == b` で一致判定）。
- 既定の電解質条件・濃度はコード内で設定されています（`TmAlign.py` 内参照）。

例
- 高感度探索（Tm しきい値をゼロ、同一率 55%以上、8 並列）
  - `python3 TmAlign.py --fasta your.fa --query ACGT... --k 3 --min-tm 0 --min-identity 40`
