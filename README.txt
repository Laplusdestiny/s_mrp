最終更新：2017/02/16 (by Toru Sumi)

### Github wikiにて更新状況を変更中
https://github.com/Laplusdestiny/s_mrp/wiki

# プログラムの使い方
## プログラムファイル
|name|中身|
|:-:|:-:|
|encmrp.c|encoderファイル|
|decmrp.c|decoderファイル|
|common.c|関数をまとめてあるファイル|
|log.c|LOGを書き出す関数をまとめてあるファイル|
|rc.c|Rangecoder本体|
|mrp.h|Headerファイル(オプション，変数系はここ)|
|run.sh|コマンドをまとめて処理してくれるシェルファイル|
|convert.sh|出力pgm,ppmをpngに変換する(run.shに一部処理が入ってる)|
|makefile|コンパイル用，コンパイルオプションなどはここで設定|
|m|コンパイル→実行を一括で(shell)|
|init_directory.sh|フォルダを再構築|
|README.txt|これ|

## オプション周り
mrp.hの中身について
### OPTIMIZE
|定義名|意味|
|:-:|:-:|
|OPT_SIDE_INFO|付加情報も含めた最適化|
|RENEW_ADC|予測器自動削除の個数を制限(墨氏)|
|PAST_ADC|従来の予測器自動削除(柴崎氏)|
|MAX_DEL_CLASS|削除個数の制限数|
|EXTRA_AUTO_DEL|自動削除に失敗してループを抜ける数|

### MULT PEAK
|定義名|意味|
|:-:|:-:|
|MULT_PEAK_MODE|多峰性確率モデルによる符号化|
|OPTIMIZE_MASK|マスクサイズの最適化|
|WIN_BSIZE|マスクサイズの切り替え単位のブロックサイズ|
|INIT_MASK|最初のマスクサイズ(基本的には1*1)|
|MAX_PEAK_NUM|多峰性確率モデルにおけるピークの最大数|

### MULT PEAK
|定義名|意味|
|:-:|:-:|
|MULT_PEAK_MODE|多峰性確率モデルによる符号化|
|OPTIMIZE_MASK|マスクサイズの最適化|
|WIN_BSIZE|マスクサイズの切り替え単位のブロックサイズ|
|INIT_MASK|最初のマスクサイズ(基本的には1*1)|
|MAX_PEAK_NUM|多峰性確率モデルにおけるピークの最大数|

### PMODEL
|定義名|意味|
|:-:|:-:|
|PM_ACCURACY|確率モデルの精度|
|CONTEXT_COST_MOUNT|特徴量を符号量で定義するか|
|CONTEXT_ERROR|特徴量を予測誤差で定義するか|
|MAX_COST_WEIGHT|符号量の重み(最大25.5で8bitでの等長符号化)|

### TEMPLATE MATCHING
|定義名|意味|
|:-:|:-:|
|TEMPLATE_MATCHING_ON|テンプレートマッチングを加えるかどうか|
|ZNCC|マッチングコストの関数としてZNCCを使用|
|MANHATTAN_SORT|事例の情報を市街地距離が近い順に更にソート|
|TEMPLATE_LOG_OUTPUT|事例の情報を吐き出す(デバックなどで時間短縮できる，0にすることで出力ファイルから事例の情報を取得)
|AREA|テンプレートの画素数|
|Y_SIZE,X_SIZE|探索範囲(x:x-X_SIZE~x+X_SIZE,y:y-Y_SIZE)|
|NAS_ACCURACY|マッチングコストの保存する精度(0.001精度)|
|MAX_DATA_SAVE|保存する事例の数|
|MAX_DATA_SAVE_DOUBLE|配列の都合上，4倍の値|
|W_GR|片側ラプラス関数の初期分散|
|W_CN|ラプラス関数の形状パラメータ|
|TEMPLATE_CLASS_NUM|テンプレートマッチングで用いる事例の最大数|
|ZSAD|マッチングコストのコスト関数としてZSADを使用(ZSADもZNCCもオフにするとSSD)|
|TEMPLATE_FLAG|予測係数として最大値を格納(予測器の先頭の係数にこれを保存している)|

### DEBUG
|定義名|意味|
|:-:|:-:|
|CHECK_TM|テンプレートマッチング周りの情報を表示|
|CHECK_TM_DETAIL|上記の詳細情報|
|CHECK_TM_WEIGHT|temp_mask_parameterでの表示|
|CHECK_DEBUG|随所で使用(エンコーダ側とデコーダ側で輝度を表示できる|
|CHECK_PMODEL|check_y,check_xと併用して特定の位置での確率モデルを表示|
|CHECK_CLASS|クラスを表示|
|CHECK_PREDICTOR|予測係数を表示|
|check_y,check_x|表示したい画素位置|
|F_NUM|一部表示において，特定の関数内での表示を行う|
|NUM_THREADS|openmp,opencv,intelコンパイラなどでの並列処理時に使用する最大スレッド数|