 ## Description


Sample code for obtaining eigenvalues using Sakurai-Sugiura (SS) method. 

## Compile 
とりあえず、sekirei (issp system B)
で動くようにしてあります。
それ以外のシステムを使うときは、
lapackのオプションと
komegaの
libとincludeを適当に変更してください。

```
sh ./com.sh
```
で実行体
SSkomega
ができます。

## Calculation
MatrixMarket方式のハミルトニアン
zvo_Ham.dat
を用意します。

./SSKomega

でそのハミルトニアンを用いてSS法を行って
指定した範囲内での、
固有値・固有ベクトルを計算します。

## SS法のパラメーター
ser parametersのところで、
コードの中(SSKomega.c)で直接指定してしまっています。

Komegaのパラメータは以下の２つです。
```
itermax: komegaのitetationの最大値
threshold : komegaの収束判定
```

周回積分は
zj=γ+ρ*exp[(2πI/nz)*(j+1/2)]
で行うようにしています。
```
gamma: 周回積分の原点
rho: 周回積分の長さ
nz: 周回積分を行う分点数
```

以下は、SS法のパラメータです。
```
nr: SS法で用いるベクトルの数
nk: SS法で計算する次数 
```

Author: Takahiro Misawa (ISSP, Univ. of Tokyo), Date: 2019/12/30
