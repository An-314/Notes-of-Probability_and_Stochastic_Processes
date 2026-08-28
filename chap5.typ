#import "@preview/scripst:1.1.2": *

= 布朗运动与极限定理

== 特征函数与高斯分布

=== 特征函数

#definition(subname: [特征函数])[
  对任意实随机变量 $X$，
  $
    phi_X(theta)=E(e^(i theta X))
    =E(cos(theta X))+i E(sin(theta X)), quad theta in RR
  $
  称为 $X$ 的特征函数。
]

由于 $abs(e^(i theta X))=1$，特征函数对一切随机变量都存在。离散与连续情形分别为
$
  phi_X(theta)=sum_k p_k e^(i theta x_k),
$
$
  phi_X(theta)=integral_(-oo)^oo e^(i theta x)f_X(x) dif x,
$
后者就是密度的 Fourier 变换。

#theorem(subname: [基本性质])[
  特征函数满足：

  - $phi_X(0)=1$，$abs(phi_X(theta))<=1$；
  - $phi_X$ 在 $RR$ 上一致连续；
  - $overline(phi_X(theta))=phi_X(-theta)$；
  - $phi_(a+b X)(theta)=e^(i a theta)phi_X(b theta)$；
  - 若 $X,Y$ 独立，则 $phi_(X+Y)(theta)=phi_X(theta)phi_Y(theta)$。
]

最后一条把独立和的卷积化为乘法，是特征函数方法的核心。

若 $E(abs(X)^k)<oo$，则可在零点求导：
$
  phi_X^(k)(0)=i^k E(X^k).
$
特别地，当二阶矩存在时，$theta arrow.r.long 0$ 有
$
  phi_X(theta)=1+i theta E(X)-theta^2 E(X^2)/2+o(theta^2).
$

常见分布的特征函数为
$
  X tilde "B"(n,p): quad phi_X(theta)=(q+p e^(i theta))^n,
$
$
  X tilde "Po"(lambda): quad phi_X(theta)=exp(lambda(e^(i theta)-1)),
$
$
  X tilde N(mu,sigma^2): quad phi_X(theta)=exp(i mu theta-sigma^2 theta^2/2).
$

#theorem(subname: [唯一性与连续性定理])[
  一个分布由其特征函数唯一决定。若 $X_n$ 的特征函数逐点收敛到在 $0$ 处连续的函数 $phi$，且 $phi$ 是随机变量 $X$ 的特征函数，则 $X_n$ 依分布收敛到 $X$；反之亦然。
]

这里“$X_n$ 依分布收敛到 $X$”是指在 $F_X$ 的所有连续点 $x$ 上，
$F_(X_n)(x) arrow.r.long F_X(x)$。

#example(subname: [泊松分布的正态近似])[
  设 $X_lambda tilde "Po"(lambda)$，令
  $Z_lambda=(X_lambda-lambda)/sqrt(lambda)$。则
  $
    phi_(Z_lambda)(theta)
    =exp(lambda(e^(i theta/sqrt(lambda))-1)-i theta sqrt(lambda)).
  $
  展开指数得其极限为 $e^(-theta^2/2)$，故
  $Z_lambda$ 依分布收敛到 $N(0,1)$。
]

=== 多维特征函数与高斯向量

对 $n$ 维随机向量 $X=(X_1,dots,X_n)^T$，定义
$
  phi_X(theta)=E(e^(i theta^T X)), quad theta in RR^n.
$
若 $Y=A X+b$，则
$
  phi_Y(theta)=e^(i theta^T b)phi_X(A^T theta).
$
向量各分量相互独立，当且仅当
$
  phi_X(theta_1,dots,theta_n)=product_(k=1)^n phi_(X_k)(theta_k).
$

#definition(subname: [多维高斯分布])[
  若 $Z=(Z_1,dots,Z_n)^T$ 的各分量独立且服从 $N(0,1)$，令
  $
    X=mu+A Z,
  $
  则称 $X$ 为高斯随机向量。其均值为 $mu$，协方差矩阵为
  $Sigma=A A^T$，记作 $X tilde N(mu,Sigma)$。
]

协方差矩阵的元素为 $Sigma_(j k)="Cov"(X_j,X_k)$，它是对称半正定矩阵。高斯向量的特征函数为
$
  phi_X(theta)=exp(i theta^T mu-1/2 theta^T Sigma theta).
$
这也可以作为高斯分布的刻画：向量 $X$ 是高斯的，当且仅当任意线性组合 $a^T X$ 都是一维正态变量。

#theorem(subname: [高斯向量中的独立性])[
  高斯随机向量的若干分量相互独立，当且仅当它们两两不相关。特别地，若协方差矩阵是对角阵，则各分量独立。
]

若 $Sigma$ 正定，$m$ 维高斯向量具有密度
$
  f_X(x)=1/((2 pi)^(m/2)sqrt(det Sigma))
  exp(-1/2 (x-mu)^T Sigma^(-1)(x-mu)).
$

== 布朗运动

=== 定义与有限维分布

#definition(subname: [标准布朗运动])[
  实值随机过程 $B=\{B_t:t>=0\}$ 若满足：

  - $B_0=0$；
  - 不相交时间区间上的增量相互独立；
  - $B_(t+s)-B_s tilde N(0,t)$；
  - 几乎所有样本轨道 $t mapsto B_t$ 连续，

  则称为标准布朗运动或 Wiener 过程。
]

漂移为 $mu$、波动率为 $sigma>0$ 的布朗运动写作
$
  X_t=x+mu t+sigma B_t.
$
它的增量满足 $X_(t+s)-X_s tilde N(mu t,sigma^2 t)$。

由于 $B_t=sum_(j=1)^k (B_(t_j)-B_(t_(j-1)))$，定义完全决定所有有限维分布。对 $0<=s<=t$，
$
  E(B_t)=0, quad D(B_t)=t,
$
$
  "Cov"(B_s,B_t)=E(B_s(B_s+B_t-B_s))=s=min(s,t).
$
因此 $(B_(t_1),dots,B_(t_n))$ 是均值为零、协方差矩阵
$Sigma_(j k)=min(t_j,t_k)$ 的高斯向量。布朗运动也是协方差函数为 $min(s,t)$ 的高斯过程。

=== Markov 性、平移与缩放

#theorem(subname: [布朗运动的基本不变性])[
  对任意 $s>=0$、$c>0$，下列过程仍是标准布朗运动：
  $
    \{B_(s+t)-B_s:t>=0\},
  $
  $
    \{c^(-1/2)B_(c t):t>=0\},
  $
  以及在 $t>0$ 上的时间反演过程 $\{t B_(1/t):t>0\}$（补上 $t=0$ 的连续值）。
]

第一条说明从任意时刻重新观察，未来增量与过去独立且规律不变，即 Markov 性的直观来源；第二条说明时间放大 $c$ 倍时，空间尺度只放大 $sqrt(c)$ 倍。

布朗运动的转移密度为
$
  p_t(x,y)=1/sqrt(2 pi t) exp(-(y-x)^2/(2t)).
$
Chapman--Kolmogorov 方程
$
  p_(s+t)(x,z)=integral_(-oo)^oo p_s(x,y)p_t(y,z) dif y
$
正是正态分布卷积的可加性。

=== 路径性质与二次变差

尽管样本轨道连续，布朗运动几乎处处不可微。形式上，在长度 $Delta t$ 的小区间中，增量量级为 $sqrt(Delta t)$，所以差商量级为 $1/sqrt(Delta t)$，不会稳定到有限导数。

对分割 $Pi:0=t_0<t_1<dots<t_n=t$，二次变差为
$
  Q_Pi=sum_(k=1)^n (B_(t_k)-B_(t_(k-1)))^2.
$
当网格宽度趋于零时，$Q_Pi$ 依概率收敛到 $t$。这常简写为
$
  (dif B_t)^2=dif t,
$
是 Itô 微积分与普通微积分不同的根源。

#example(subname: [首次到达的反射原理])[
  令 $M_t=max_(0<=s<=t)B_s$。对 $a>0$，反射首次到达 $a$ 之后的轨道，得到
  $
    P(M_t>=a)=2P(B_t>=a)=2(1-Phi(a/sqrt(t))).
  $
  因而首次到达时刻 $tau_a=inf\{t:B_t=a\}$ 满足
  $P(tau_a<=t)=2(1-Phi(a/sqrt(t)))$。
]

== 概率收敛与大数定律

=== 依概率收敛

#definition(subname: [依概率收敛])[
  若对每个 $epsilon>0$，
  $
    P(abs(X_n-X)>epsilon) arrow.r.long 0,
  $
  则称 $X_n$ 依概率收敛到 $X$，记作 $X_n arrow.r.long^P X$。
]

依概率收敛推出依分布收敛，反之一般不成立。若极限是常数，则依分布收敛与依概率收敛等价。

#theorem(subname: [Chebyshev 大数定律])[
  若 $X_1,X_2,dots$ 相互独立，$E(X_k)=mu_k$，且方差一致有界 $D(X_k)<=C$，则
  $
    1/n sum_(k=1)^n (X_k-mu_k) arrow.r.long^P 0.
  $
]

_证明._ 独立性给出
$
  D(1/n sum_(k=1)^n X_k)=1/n^2 sum_(k=1)^n D(X_k)<=C/n.
$
应用 Chebyshev 不等式即可。

独立同分布且 $E(X_1)=mu$、$D(X_1)<oo$ 时，得到弱大数定律
$
  overline(X)_n=1/n sum_(k=1)^n X_k arrow.r.long^P mu.
$

#theorem(subname: [Bernoulli 大数定律])[
  在独立重复试验中，事件 $A$ 的概率为 $p$，令 $N_n$ 为前 $n$ 次中 $A$ 的发生次数，则
  $
    N_n/n arrow.r.long^P p.
  $
  所以频率在大量试验中稳定到概率。
]

更一般的 Markov 条件是
$
  1/n^2 D(sum_(k=1)^n X_k) arrow.r.long 0,
$
它已足以保证样本均值减去均值的平均依概率趋于零，不必要求独立。

#example(subname: [Monte Carlo 积分])[
  若 $U_k$ 独立同服从 $"U"(0,1)$，且 $g$ 可积，则
  $
    1/n sum_(k=1)^n g(U_k) arrow.r.long^P integral_0^1 g(x) dif x.
  $
  因而可通过均匀随机点上的函数平均估计积分。若方差有限，均方误差按 $1/n$ 衰减，典型误差按 $1/sqrt(n)$ 衰减。
]

== 中心极限定理

#theorem(subname: [Lindeberg--Lévy 中心极限定理])[
  设 $X_1,X_2,dots$ 独立同分布，$E(X_k)=mu$，$D(X_k)=sigma^2 in (0,oo)$。则
  $
    S_n=(sum_(k=1)^n X_k-n mu)/(sigma sqrt(n))
  $
  依分布收敛到 $N(0,1)$，即对任意 $x$，
  $
    lim_(n arrow.r.long oo)P(S_n<=x)=Phi(x).
  $
]

_特征函数证明._ 令 $Y_k=(X_k-mu)/sigma$，则
$
  phi_Y(theta)=1-theta^2/2+o(theta^2).
$
独立性给出
$
  phi_(S_n)(theta)=(phi_Y(theta/sqrt(n)))^n
  =(1-theta^2/(2 n)+o(1/n))^n arrow.r.long e^(-theta^2/2).
$
由特征函数连续性定理，极限分布为标准正态。

#theorem(subname: [棣莫弗--Laplace 定理])[
  若 $X_n tilde "B"(n,p)$，$q=1-p$，则
  $
    (X_n-n p)/sqrt(n p q) arrow.r.long^d N(0,1).
  $
  因而大样本下
  $
    P(a<=X_n<=b) approx
    Phi((b+1/2-n p)/sqrt(n p q))
    -Phi((a-1/2-n p)/sqrt(n p q)),
  $
  其中 $plus.minus 1/2$ 是把离散区间近似为连续区间时的连续性校正。
]

大数定律回答“样本平均是否靠近总体均值”，中心极限定理进一步回答“偏差按 $sqrt(n)$ 放大后呈现什么形状”：前者给稳定性，后者给近似分布和误差尺度。

#example(subname: [独立负载的正态近似])[
  若 $n$ 台设备独立开机，单台开机概率 $p$、开机功率 $w$，总功率为 $W=w X$，其中 $X tilde "B"(n,p)$。于是
  $
    E(W)=w n p, quad D(W)=w^2 n p q.
  $
  当 $n p$ 与 $n q$ 都足够大时，可以用相应正态分布估计超出供电容量的概率，并使用连续性校正。
]

#example(subname: [舍入误差的累积])[
  若各次独立舍入误差 $E_k tilde "U"(-1/2,1/2)$，则
  $E(E_k)=0$、$D(E_k)=1/12$。$n$ 次误差之和近似服从
  $
    N(0,n/12),
  $
  说明独立误差的标准差按 $sqrt(n)$ 而不是 $n$ 增长。
]
