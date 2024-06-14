このプログラムではc++向けの行列演算ライブラリであるEigenを使用しています。
まず最初にEigenをインストールします。
```
sudo apt install libeigen3-dev
```

次に以下のコマンドでコンパイルしてください。
```
g++ -I /usr/include/eigen3/ main.cpp
```

実行結果は以下の通りです。
```
$ ./a.out
Jacobi Method in EQ1
tol: 1e-06
Epoch: 45
1
2
3

Gauss-Seidel Method in EQ1
tol: 1e-06
Epoch: 21
1
2
3

Jacobi Method in EQ2
tol: 1e-06
Epoch: 61
-5
2
4
1
-3

Gauss-Seidel Method in EQ2
tol: 1e-06
Epoch: 12
-5
2
4
1
-3
```