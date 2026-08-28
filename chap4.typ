#import "@preview/scripst:1.1.2": *

= 连续型随机变量

== 分布函数与概率密度

#definition(subname: [分布函数])[
  对任意实随机变量 $X$，
  $
    F_X(x)=P(X<=x), quad x in RR
  $
  称为 $X$ 的分布函数。
]

分布函数同时适用于离散与连续随机变量，并满足：

- $0<=F_X(x)<=1$，且单调不减；
- $lim_(x arrow.r.long -oo)F_X(x)=0$，$lim_(x arrow.r.long +oo)F_X(x)=1$；
- $F_X$ 右连续；
- $P(a<X<=b)=F_X(b)-F_X(a)$；
- $P(X=x)=F_X(x)-F_X(x^-)$。

因此，分布函数的跳跃恰好对应点质量。连续型随机变量在每个单点上的概率为零，但区间事件仍可具有正概率；“概率为零”不等于事件不可能发生。

#definition(subname: [连续型随机变量与密度])[
  若存在非负可积函数 $f_X$，使
  $
    F_X(x)=integral_(-oo)^x f_X(t) dif t,
  $
  则称 $X$ 为连续型随机变量，$f_X$ 为概率密度。
]

密度满足
$
  f_X(x)>=0, quad integral_(-oo)^oo f_X(x) dif x=1,
$
$
  P(X in A)=integral_A f_X(x) dif x.
$
在 $f_X$ 连续处有 $F_X'(x)=f_X(x)$。密度值可以大于 $1$；真正的概率是曲线下的面积。

#example(subname: [由常数确定密度])[
  若 $f(x)=c x(1-x) I_((0,1))(x)$，归一化给出
  $
    1=c integral_0^1 x(1-x) dif x=c/6,
  $
  故 $c=6$。于是对 $0<=a<b<=1$，
  $P(a<X<=b)=integral_a^b 6x(1-x) dif x$。
]

== 数学期望、方差与不等式

若积分绝对收敛，则
$
  E(X)=integral_(-oo)^oo x f_X(x) dif x.
$
更一般地，对可测函数 $g$，
$
  E(g(X))=integral_(-oo)^oo g(x)f_X(x) dif x.
$
这称为无意识统计学家定律：求函数的期望时无需先求 $g(X)$ 的分布。

方差为
$
  D(X)=E((X-E(X))^2)=E(X^2)-E(X)^2.
$
期望的线性性不要求独立；而独立随机变量满足
$
  E(X Y)=E(X)E(Y), quad D(X+Y)=D(X)+D(Y).
$

#theorem(subname: [Cauchy--Schwarz 不等式])[
  若二阶矩存在，则
  $
    abs(E(X Y))^2<=E(X^2)E(Y^2).
  $
  取 $Y=1$ 可得 $abs(E(X))<=sqrt(E(X^2))$。
]

#theorem(subname: [Markov 与 Chebyshev 不等式])[
  对非负随机变量 $Y$ 和 $a>0$，
  $
    P(Y>=a)<=E(Y)/a.
  $
  令 $Y=(X-E(X))^2$，便得到
  $
    P(abs(X-E(X))>=epsilon)<=D(X)/epsilon^2.
  $
]

== 常见连续分布

=== 均匀分布

若
$
  f_X(x)=cases(1/(b-a), & a<=x<=b, 0, & "其他"),
$
则记 $X tilde "U"(a,b)$。其分布函数在 $[a,b]$ 上线性增长，并且
$
  E(X)=(a+b)/2, quad D(X)=(b-a)^2/12.
$
若 $U tilde "U"(0,1)$，则 $a+(b-a)U tilde "U"(a,b)$。

=== 指数分布与 Gamma 分布

若
$
  f_X(x)=lambda e^(-lambda x)I_((0,+oo))(x), quad lambda>0,
$
则记 $X tilde "Exp"(lambda)$。它的分布函数与生存函数为
$
  F_X(x)=1-e^(-lambda x), quad P(X>x)=e^(-lambda x) quad (x>=0),
$
且
$
  E(X)=1/lambda, quad D(X)=1/lambda^2.
$

#theorem(subname: [无记忆性])[
  若 $X tilde "Exp"(lambda)$，则对 $s,t>=0$，
  $
    P(X>s+t | X>s)=P(X>t).
  $
  反过来，在适当正则性条件下，非负连续分布中只有指数分布具有无记忆性。
]

Gamma 函数定义为
$
  Gamma(alpha)=integral_0^oo u^(alpha-1)e^(-u) dif u,
$
满足 $Gamma(alpha+1)=alpha Gamma(alpha)$、$Gamma(n)=(n-1)!$、$Gamma(1/2)=sqrt(pi)$。若
$
  f(x)=lambda^alpha/Gamma(alpha) x^(alpha-1)e^(-lambda x)I_((0,+oo))(x),
$
则称 $X$ 服从形状 $alpha$、率 $lambda$ 的 Gamma 分布，并有
$
  E(X)=alpha/lambda, quad D(X)=alpha/lambda^2.
$
整数形状时它是泊松过程中第 $alpha$ 次到达的等待时间。

=== 正态分布

若
$
  f_X(x)=1/(sqrt(2 pi)sigma) exp(-(x-mu)^2/(2sigma^2)), quad sigma>0,
$
则记 $X tilde N(mu,sigma^2)$。通过二维极坐标可证明高斯积分
$
  integral_(-oo)^oo e^(-x^2/2) dif x=sqrt(2 pi),
$
从而密度归一化。并且
$
  E(X)=mu, quad D(X)=sigma^2.
$

标准正态密度与分布函数记作
$
  phi(x)=1/sqrt(2 pi)e^(-x^2/2), quad Phi(x)=integral_(-oo)^x phi(t) dif t.
$
有 $phi(-x)=phi(x)$、$Phi(-x)=1-Phi(x)$。标准化
$
  Z=(X-mu)/sigma tilde N(0,1)
$
给出
$
  P(a<=X<=b)=Phi((b-mu)/sigma)-Phi((a-mu)/sigma).
$
常用的经验数值为
$
  P(abs(X-mu)<=sigma) approx 0.6827,
$
$
  P(abs(X-mu)<=2sigma) approx 0.9545, quad
  P(abs(X-mu)<=3sigma) approx 0.9973.
$

== 二维连续型随机变量

=== 联合分布、边缘分布

对随机向量 $(X,Y)$，联合分布函数为
$
  F_(X,Y)(x,y)=P(X<=x,Y<=y).
$
若存在联合密度 $f_(X,Y)$，则
$
  F_(X,Y)(x,y)=integral_(-oo)^x integral_(-oo)^y f_(X,Y)(u,v) dif v dif u.
$
对平面区域 $D$，
$
  P((X,Y) in D)=integral.double_D f_(X,Y)(x,y) dif x dif y.
$
边缘密度由另一变量积分消去：
$
  f_X(x)=integral_(-oo)^oo f_(X,Y)(x,y) dif y,
$
$
  f_Y(y)=integral_(-oo)^oo f_(X,Y)(x,y) dif x.
$

#example(subname: [区域上的均匀分布])[
  若 $D subset RR^2$ 面积为 $|D|$，区域上的均匀分布密度为
  $f(x,y)=I_D(x,y)/|D|$。对单位圆盘 $x^2+y^2<=1$，密度为 $1/pi$；其边缘密度
  $
    f_X(x)=cases(2sqrt(1-x^2)/pi, & abs(x)<=1, 0, & "其他")
  $
  并不是区间上的均匀分布。
]

=== 二维正态分布

二维正态分布 $N(mu_1,mu_2,sigma_1^2,sigma_2^2,rho)$ 的密度为
$
  f(x,y)=1/(2 pi sigma_1 sigma_2 sqrt(1-rho^2))
  exp(-1/(2(1-rho^2))
    ((x-mu_1)^2/sigma_1^2-2rho (x-mu_1)(y-mu_2)/(sigma_1 sigma_2)
    +(y-mu_2)^2/sigma_2^2)),
$
其中 $abs(rho)<1$。边缘分布分别为
$
  X tilde N(mu_1,sigma_1^2), quad Y tilde N(mu_2,sigma_2^2).
$

=== 独立性、协方差和相关系数

$X,Y$ 相互独立，当且仅当对所有 $x,y$，
$
  F_(X,Y)(x,y)=F_X(x)F_Y(y).
$
若有密度，则等价于几乎处处
$
  f_(X,Y)(x,y)=f_X(x)f_Y(y).
$
独立时，对适当函数 $g,h$，$E(g(X)h(Y))=E(g(X))E(h(Y))$。

协方差与相关系数为
$
  "Cov"(X,Y)=E((X-E(X))(Y-E(Y)))=E(X Y)-E(X)E(Y),
$
$
  rho_(X Y)="Cov"(X,Y)/(sqrt(D(X))sqrt(D(Y))).
$
Cauchy--Schwarz 不等式保证 $abs(rho_(X Y))<=1$。$abs(rho)=1$ 当且仅当 $Y=a X+b$ 几乎必然成立；$rho=0$ 只表示线性不相关，一般不推出独立。但对联合正态变量，零协方差等价于独立。

#example(subname: [不相关但不独立])[
  若 $(X,Y)$ 在单位圆盘上均匀分布，对称性给出
  $E(X)=E(Y)=E(X Y)=0$，所以相关系数为零；然而联合密度的支撑是圆盘而不是两个边缘支撑的直积，故二者不独立。
]

== 条件分布与条件期望

=== 条件密度

当 $f_Y(y)>0$ 时，定义
$
  f_(X|Y)(x|y)=f_(X,Y)(x,y)/f_Y(y).
$
它对 $x$ 是归一化密度。相应地
$
  F_(X|Y)(x|y)=integral_(-oo)^x f_(X|Y)(u|y) dif u.
$
这可以从窄条事件 $y<Y<=y+h$ 的条件概率令 $h arrow.r.long 0+$ 得到。

对于二维正态分布，配方可得
$
  X|Y=y tilde N(mu_1+rho sigma_1/sigma_2 (y-mu_2),
    sigma_1^2(1-rho^2)).
$

=== 条件期望

条件期望是 $Y$ 的函数：
$
  E(X|Y=y)=integral_(-oo)^oo x f_(X|Y)(x|y) dif x,
$
$
  E(X|Y)=g(Y), quad g(y)=E(X|Y=y).
$
更一般地，
$
  E(h(X)|Y=y)=integral h(x)f_(X|Y)(x|y) dif x.
$

#theorem(subname: [全期望与全方差])[
  若相应矩存在，则
  $
    E(X)=E(E(X|Y)),
  $
  $
    D(X)=E(D(X|Y))+D(E(X|Y)).
  $
  后一式把总波动分解为组内波动与组间波动。
]

条件期望还满足线性性、单调性以及“已知量可以提出”：
$
  E(a X+b Z|Y)=a E(X|Y)+b E(Z|Y),
$
$
  E(h(Y)X|Y)=h(Y)E(X|Y).
$
若 $X$ 与 $Y$ 独立，则 $E(X|Y)=E(X)$。

#theorem(subname: [最佳均方预测])[
  在所有形如 $phi(Y)$ 的平方可积估计中，$E(X|Y)$ 使
  $E((X-phi(Y))^2)$ 最小。事实上
  $
    E((X-phi(Y))^2)
    =E((X-E(X|Y))^2)+E((E(X|Y)-phi(Y))^2).
  $
  等价地，残差 $X-E(X|Y)$ 与一切 $Y$ 的平方可积函数正交。
]

== 随机变量函数的分布

=== 一维变换

若 $Y=g(X)$ 且 $g$ 严格单调、可微，反函数为 $x=g^(-1)(y)$，则
$
  f_Y(y)=f_X(g^(-1)(y)) abs((dif)/(dif y)g^(-1)(y)).
$
若同一个 $y$ 有多个原像 $x_i$，则对各分支求和：
$
  f_Y(y)=sum_i f_X(x_i)/abs(g'(x_i)).
$

#example(subname: [线性变换与平方变换])[
  对 $Y=a X+b$（$a!=0$），
  $f_Y(y)=1/abs(a) f_X((y-b)/a)$。

  对 $Y=X^2$，$y>0$ 有两个原像 $plus.minus sqrt(y)$，故
  $
    f_Y(y)=(f_X(sqrt(y))+f_X(-sqrt(y)))/(2sqrt(y)).
  $
]

=== 多维变量变换

若 $(U,V)=T(X,Y)$ 是一一、光滑变换，逆变换为 $(x(u,v),y(u,v))$，则
$
  f_(U,V)(u,v)=f_(X,Y)(x(u,v),y(u,v))
  abs(det(partial(x,y)/partial(u,v))).
$
雅可比行列式描述面积元在坐标变换下的伸缩。

=== 随机变量之和与卷积

对任意具有联合密度的 $X,Y$，令 $S=X+Y$，则
$
  f_S(s)=integral_(-oo)^oo f_(X,Y)(x,s-x) dif x.
$
若 $X,Y$ 独立，便成为卷积公式
$
  f_S(s)=integral_(-oo)^oo f_X(x)f_Y(s-x) dif x=(f_X ast f_Y)(s).
$

#example(subname: [正态分布的可加性])[
  若 $X tilde N(mu_1,sigma_1^2)$、$Y tilde N(mu_2,sigma_2^2)$ 且独立，则
  $
    X+Y tilde N(mu_1+mu_2,sigma_1^2+sigma_2^2).
  $
  可由卷积直接计算，也可在下一章用特征函数一行得到。
]

#example(subname: [两个独立指数等待时间的先后])[
  若 $X_i tilde "Exp"(lambda_i)$ 独立，则
  $
    P(X_1<X_2)=integral_0^oo lambda_1 e^(-lambda_1 x)e^(-lambda_2 x) dif x
    =lambda_1/(lambda_1+lambda_2).
  $
  且 $min(X_1,X_2) tilde "Exp"(lambda_1+lambda_2)$。
]
