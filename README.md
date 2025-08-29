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

```
63.29	gi|15079186|ref|NC_002333.2|tax|7955|Danio	296	316	+	95.2	GTCGGTAAAACTCGTGCCAGC	GCCGGTAAAACTCGTGCCAGC
53.49	gi|15079186|ref|NC_002333.2|tax|7955|Danio	297	316	+	90.5	GTCGGTAAAACTCGTGCCAGC	-CCGGTAAAACTCGTGCCAGC
36.21	gi|15079186|ref|NC_002333.2|tax|7955|Danio	3437	3456	-	66.7	GCTGGCACGAGTTTTACCGAC	GCTGGCA-GAAACAAACCGAG
63.29	gi|11545708|ref|NC_002616.1|tax|41697|Sardinops	298	318	+	95.2	GTCGGTAAAACTCGTGCCAGC	GCCGGTAAAACTCGTGCCAGC
53.49	gi|11545708|ref|NC_002616.1|tax|41697|Sardinops	299	318	+	90.5	GTCGGTAAAACTCGTGCCAGC	-CCGGTAAAACTCGTGCCAGC
32.15	gi|11545708|ref|NC_002616.1|tax|41697|Sardinops	9732	9753	-	59.1	GCTGGCAC-GAGTTTTACCGAC	GCTCCCACAGATGAACCCCGAC
31.01	gi|11545708|ref|NC_002616.1|tax|41697|Sardinops	3890	3910	+	71.4	GTCGGTAAAACTCGTGCCAGC	GGGGGTTAAACTCCCCCCAGC
63.29	gi|12248136|ref|NC_002646.1|tax|59291|Coregonus	293	313	+	95.2	GTCGGTAAAACTCGTGCCAGC	GCCGGTAAAACTCGTGCCAGC
53.49	gi|12248136|ref|NC_002646.1|tax|59291|Coregonus	294	313	+	90.5	GTCGGTAAAACTCGTGCCAGC	-CCGGTAAAACTCGTGCCAGC
34.88	gi|12248136|ref|NC_002646.1|tax|59291|Coregonus	13601	13621	+	66.7	GTCGGTAAAACTCGTGCCAGC	CTCGGTCAGGCCATTGCCAGC
63.29	gi|12248150|ref|NC_002648.1|tax|81385|Polymixia	297	317	+	95.2	GTCGGTAAAACTCGTGCCAGC	GCCGGTAAAACTCGTGCCAGC
53.49	gi|12248150|ref|NC_002648.1|tax|81385|Polymixia	298	317	+	90.5	GTCGGTAAAACTCGTGCCAGC	-CCGGTAAAACTCGTGCCAGC
34.58	gi|12248150|ref|NC_002648.1|tax|81385|Polymixia	1319	1340	+	72.7	GTCGGTAAAACTCGTGC-CAGC	GCCAGTAAAACTCAAGCACAGT
30.40	gi|12248150|ref|NC_002648.1|tax|81385|Polymixia	10331	10351	-	57.1	GCTGGCACGAGTTTTACCGAC	GCACACACGGAACAGACCGAC
63.29	gi|12248164|ref|NC_002647.1|tax|91975|Diplophos	294	314	+	95.2	GTCGGTAAAACTCGTGCCAGC	GCCGGTAAAACTCGTGCCAGC
53.49	gi|12248164|ref|NC_002647.1|tax|91975|Diplophos	295	314	+	90.5	GTCGGTAAAACTCGTGCCAGC	-CCGGTAAAACTCGTGCCAGC
```

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

## TmBLAST — 複数クエリ + 縮重塩基対応

概要
- 複数のクエリ配列を一括で検索します（`--query` の複数指定、またはファイル入力）。
- クエリ配列に含まれる IUPAC 縮重塩基（R,Y,S,W,K,M,B,D,H,V,N）を A/C/G/T の全組み合わせに展開して TmAlign ロジックで検索します。

使い方
- 単純な例（2 クエリを指定）
  - `python3 TmBLAST.py --fasta target.fa --query ATGCRN --query TTYGGA --k 4`
- クエリファイルから（FASTA もしくは1行1クエリ、任意で`id\tseq`）
  - `python3 TmBLAST.py --fasta target.fa --queries-file queries.fa --k 4`

主要オプション
- `--fasta`: 参照 FASTA
- `--query`: クエリ配列（複数回指定可、IUPAC 可）
- `--queries-file`: クエリファイル（FASTA か、`id\tseq` または `seq` の行）
- `--k`, `--min-tm`, `--min-identity`, `--cpus`: `TmAlign.py` と同じ意味
- `--max-expansions`: 1 クエリあたりの展開上限（既定: 100000）

出力形式（タブ区切り）
- 列: `query_id  expanded_query  Tm  contig  start  end  strand  identity(%)  query_align  db_align`
- `query_id`: クエリ名（FASTA ヘッダ、または `q1`, `q2`, ...）
- `expanded_query`: 縮重展開後の A/C/G/T のみの配列

注意
- 縮重の展開は指数的に増大することがあります。既定では `--max-expansions` を超えると打ち切ります（標準エラーに警告）。
- 出力は各展開クエリごとに DB ヒット（座標・向き・整列）で重複排除しています。

例
- 付属のサンプルで実行:
  - `python3 TmBLAST.py --fasta test-mito.fa --queries-file test-primer.fa --k 4 --min-tm 0 --min-identity 40`
  - 出力例（先頭数行）:
```
MiFish-F	GTCGGTAAAACTCGTGCCAGC	63.29	gi|NC_000860.1|ref|NC_000860.1|tax|8038|Salvelinus	293	313	+	95.2	GTCGGTAAAACTCGTGCCAGC	GCCGGTAAAACTCGTGCCAGC
MiFish-F	GTCGGTAAAACTCGTGCCAGC	53.49	gi|NC_000860.1|ref|NC_000860.1|tax|8038|Salvelinus	294	313	+	90.5	GTCGGTAAAACTCGTGCCAGC	-CCGGTAAAACTCGTGCCAGC
MiFish-R	CATAGTGGGGTATCTAATCCCAGTTTG	66.49	gi|NC_000860.1|ref|NC_000860.1|tax|8038|Salvelinus	484	510	-	100.0	CAAACTGGGATTAGATACCCCACTATG	CAAACTGGGATTAGATACCCCACTATG
```
