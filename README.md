# SparseMatrix

- 自己实现的一个稀疏矩阵模板类，构造函数参数为矩阵行列数
- 定义矩阵之后需要初始化，初始化方式有三种，分别为：
    - 以```vector<pair>```或者```vector<tuple>```给出矩阵轮廓
    - 初始化成对角矩阵，默认单位矩阵
    - 从文件初始化
- 支持运算符有： +, +=, -, -=, *=, =, ==, <<
- 支持矩阵与向量乘，矩阵与稠密矩阵乘，详见函数```MultiplyVector```,```MultiplyVectorAdd```, ```TransposeMultiplyVector```, ```TransposeMultiplyVectorAdd```, ```MultiplyMatrix```, ```MultiplyMatrixAdd```
- 此外，还支持操作有： ```set(i, j)```, ```get(i, j)```, ```CheckEntry```, ```CheckSysmetric()```, ```ResetZero()```, ```Print```等，详细见代码。