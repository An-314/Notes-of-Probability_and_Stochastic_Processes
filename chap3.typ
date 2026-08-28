#import "@preview/scripst:1.1.2": *

= 随机徘徊与泊松过程

本章从离散时间的随机徘徊出发，引入随机过程、状态空间、样本轨道和有限维分布；随后以概率母函数研究泊松分布，并建立连续时间的泊松过程。两类过程都是“独立平稳增量”的基本模型。

== 随机徘徊

=== 一维简单随机徘徊

不断独立投掷一枚硬币，正面概率为 $p$，反面概率为 $q=1-p$。令
$
  Y_n = cases(1, "第" n "次为正面", 0, "第" n "次为反面"),
  quad Z_n=2Y_n-1,
$
则 $P(Z_n=1)=p$、$P(Z_n=-1)=q$。独立同分布序列 $\{Z_n:n>=1\}$ 称为 Bernoulli 序列。若粒子从 $X_0=x$ 出发，每次按 $Z_n$ 移动一步，则
$
  X_n=X_0+sum_(k=1)^n Z_k
$
称为一维简单随机徘徊；$p=q=1/2$ 时称为对称简单随机徘徊。

对固定的一串试验结果 $omega$，函数 $n mapsto X_n(omega)$ 称为一条样本轨道；所有可能位置构成状态空间 $ZZ$。

#theorem(subname: [位置的分布])[
  设 $X_0=0$。在 $n$ 步后到达 $s$，等价于向右走了 $(n+s)/2$ 步。因此
  $
    P(X_n=s)=cases(
      binom(n,(n+s)/2)p^((n+s)/2)q^((n-s)/2), & n>=abs(s) "且" n+s "为偶数", \
      0, & "其他情形".
    )
  $
  特别地，若 $S_n=(X_n+n)/2$，则 $S_n tilde "B"(n,p)$。
]

由于不相交时段的步长和相互独立，对 $m>n$ 有
$
  P(X_n=r,X_m=s)=P(X_n=r)P(X_(m-n)=s-r).
$
更一般地，$0=n_0<n_1<dots<n_k$ 时，
$
  P(X_(n_1)=i_1,dots,X_(n_k)=i_k)
  =product_(j=1)^k P(X_(n_j-n_(j-1))=i_j-i_(j-1)),
$
其中 $i_0=X_0$。这给出了随机徘徊的全部有限维分布。

由 $E(Z_k)=p-q=2p-1$、$D(Z_k)=4 p q$，
$
  E(X_n)=x+n(2p-1), quad D(X_n)=4 n p q.
$
当 $m>=n$ 时，$X_m=X_n+sum_(k=n+1)^m Z_k$，后一个增量与 $X_n$ 独立，故
$
  "Cov"(X_n,X_m)=D(X_n)=4 n p q.
$
于是一般地
$
  "Cov"(X_n,X_m)=4 p q min(n,m).
$

=== 首次到达与吸收边界

令
$
  T_a=inf\{n>=0:X_n=a\}
$
为首次到达 $a$ 的时刻；若永不到达，则约定 $T_a=+oo$。首次返回原点的时刻通常写作 $T_0^+=inf\{n>=1:X_n=0\}$。

#example(subname: [赌徒破产或带吸收壁的徘徊])[
  粒子从 $i in \{0,1,dots,N\}$ 出发，在 $0,N$ 处吸收。记先到 $N$ 的概率为
  $u_i=P_i(T_N<T_0)$。一步分析给出
  $
    u_i=p u_(i+1)+q u_(i-1), quad u_0=0, quad u_N=1.
  $
  解差分方程得
  $
    u_i=cases(
      i/N, & p=q=1/2, \
      (1-(q/p)^i)/(1-(q/p)^N), & p!=q.
    )
  $
  因而先到 $0$ 的概率为 $1-u_i$。对称情形中，吸收前平均步数 $e_i$ 满足
  $e_i=1+(e_(i-1)+e_(i+1))/2$，边界 $e_0=e_N=0$，故 $e_i=i(N-i)$。
]

== 随机过程的基本语言

#definition(subname: [随机过程])[
  在概率空间 $(Omega,cal(F),P)$ 上，一族以 $t in T$ 为指标的随机变量
  $X=\{X_t:t in T\}$ 称为随机过程。$T$ 是指标集，通常解释为时间；固定 $omega$ 后的函数 $t mapsto X_t(omega)$ 是一条样本轨道。
]

若 $T$ 可数，则称离散参数过程；若 $T$ 是实区间，则称连续参数过程。过程的取值集合称为状态空间，因此可以按时间和状态分别离散或连续：随机徘徊是离散时间、离散状态；泊松过程是连续时间、离散状态；布朗运动是连续时间、连续状态。

#definition(subname: [独立增量与平稳增量])[
  若任意互不相交区间上的增量相互独立，则称过程具有独立增量。若 $X_(t+h)-X_t$ 的分布只依赖区间长度 $h$ 而不依赖起点 $t$，则称具有平稳增量。
]

随机徘徊满足独立平稳增量。二维随机徘徊也可用独立步长向量 $Z_n$ 定义，例如每一步等概率取 $(1,0),(-1,0),(0,1),(0,-1)$。

== 泊松分布与概率母函数

=== 泊松分布的回顾

若
$
  P(X=k)=e^(-lambda)lambda^k/k!, quad k=0,1,2,dots,
$
则记 $X tilde "Po"(lambda)$，并有 $E(X)=D(X)=lambda$。

#theorem(subname: [泊松极限定理])[
  若 $n arrow.r.long +oo$、$p_n arrow.r.long 0$ 且 $n p_n arrow.r.long lambda>0$，则对固定 $k$，
  $
    binom(n,k)p_n^k(1-p_n)^(n-k) arrow.r.long e^(-lambda)lambda^k/k!.
  $
  因此“大量独立机会、每次小概率发生”的计数常可近似为泊松分布。
]

#theorem(subname: [可加性与稀疏化])[
  若独立随机变量 $X_i tilde "Po"(lambda_i)$，则
  $sum_i X_i tilde "Po"(sum_i lambda_i)$。若对 $X tilde "Po"(lambda)$ 个对象各自独立地以概率 $r$ 保留，保留数 $Y$ 满足 $Y tilde "Po"(lambda r)$；被舍弃数与保留数还相互独立。
]

可加性可由卷积直接验证；稀疏化可在给定 $X=n$ 后使用二项分布，再对 $n$ 求和。

=== 概率母函数

#definition(subname: [概率母函数])[
  对非负整数值随机变量 $X$，定义
  $
    G_X(s)=E(s^X)=sum_(k=0)^oo P(X=k)s^k.
  $
]

在收敛范围内，
$
  P(X=k)=G_X^(k)(0)/k!, quad E(X)=G_X'(1),
$
$
  D(X)=G_X''(1)+G_X'(1)-G_X'(1)^2.
$
若 $X_1,dots,X_n$ 相互独立，则
$
  G_(X_1+dots+X_n)(s)=product_(j=1)^n G_(X_j)(s).
$

常见母函数为
$
  G_("B"(n,p))(s)=(q+p s)^n,
$
$
  G_("Po"(lambda))(s)=exp(lambda(s-1)),
$
$
  G_("Ge"(p))(s)=(p s)/(1-q s)
$
（这里几何分布取值 $1,2,dots$）。由乘法性质立即得到二项分布和泊松分布的可加性；$n$ 个几何等待时间之和服从负二项分布
$
  P(S_n=k)=binom(k-1,n-1)p^n q^(k-n), quad k>=n.
$

#example(subname: [分枝过程中的母函数复合])[
  一个细胞产生后代数 $X$，母函数为 $g(s)=E(s^X)$。若各后代独立地按同一规则再繁殖一代，则两代后总数的母函数为
  $
    G_2(s)=E(g(s)^X)=g(g(s)).
  $
  迭代 $n$ 代得到 $G_n=g compose dots compose g$。这是“先条件化，再取期望”在母函数中的典型用法。
]

== 泊松过程

=== 定义与计数分布

#definition(subname: [泊松过程])[
  非负整数值过程 $N=\{N_t:t>=0\}$ 若满足：

  - $N_0=0$；
  - 不相交时间区间上的增量相互独立；
  - 增量平稳；
  - 当 $h arrow.r.long 0+$ 时，
    $P(N_(t+h)-N_t=1)=lambda h+o(h)$，
    $P(N_(t+h)-N_t>=2)=o(h)$，

  则称为强度 $lambda$ 的泊松过程。
]

最后一条称为普通性。由概率归一性还有
$
  P(N_(t+h)-N_t=0)=1-lambda h+o(h).
$

#theorem(subname: [泊松过程的增量分布])[
  若 $N$ 是强度为 $lambda$ 的泊松过程，则
  $
    N_(t+s)-N_t tilde "Po"(lambda s).
  $
]

_推导._ 令 $p_k(t)=P(N_t=k)$。独立增量与局部条件给出
$
  p_0(t+h)=p_0(t)(1-lambda h)+o(h),
$
所以 $p_0'(t)=-lambda p_0(t)$，即 $p_0(t)=e^(-lambda t)$。对 $k>=1$，
$
  p_k'(t)=-lambda p_k(t)+lambda p_(k-1)(t).
$
结合初值 $p_k(0)=0$ 递推得到
$
  p_k(t)=e^(-lambda t)(lambda t)^k/k!.
$
平稳增量给出一般区间的结论。反过来，“$N_0=0$、独立增量、增量服从 $"Po"(lambda t)$”也可作为等价定义。

由分布立即有
$
  E(N_t)=D(N_t)=lambda t,
$
而对 $s<=t$，写成 $N_t=N_s+(N_t-N_s)$，故
$
  "Cov"(N_s,N_t)=D(N_s)=lambda s,
$
即一般地 $"Cov"(N_s,N_t)=lambda min(s,t)$。

=== 条件分布与有限维分布

#theorem(subname: [给定总数后的均匀散布])[
  对 $0<=s<=t$，
  $
    P(N_s=k | N_t=n)=binom(n,k)(s/t)^k(1-s/t)^(n-k).
  $
  换言之，在已知 $[0,t]$ 内恰有 $n$ 次到达后，每次到达落入 $[0,s]$ 的方式如同独立概率 $s/t$ 的选择。
]

对 $0=t_0<t_1<dots<t_m$，各增量独立，因此
$
  P(N_(t_1)=n_1,dots,N_(t_m)=n_m)
  =product_(j=1)^m e^(-lambda(t_j-t_(j-1)))
   (lambda(t_j-t_(j-1)))^(n_j-n_(j-1))/(n_j-n_(j-1))!,
$
其中 $0=n_0<=n_1<=dots<=n_m$；若不满足单调性，概率为零。

=== 到达间隔与到达时刻

令 $T_n=inf\{t:N_t>=n\}$ 为第 $n$ 次到达时刻，$W_n=T_n-T_(n-1)$ 为相邻到达间隔。

#theorem(subname: [指数等待时间])[
  $W_1,W_2,dots$ 相互独立且同服从 $"Exp"(lambda)$。因此
  $
    T_n=sum_(k=1)^n W_k
  $
  服从 Gamma/Erlang 分布，密度为
  $
    f_(T_n)(t)=cases(lambda^n t^(n-1)e^(-lambda t)/(n-1)!, & t>0, 0, & t<=0),
  $
  且 $E(T_n)=n/lambda$、$D(T_n)=n/lambda^2$。
]

关键关系是 $\{T_n<=t\}=\{N_t>=n\}$。例如
$
  P(W_1>t)=P(N_t=0)=e^(-lambda t).
$

=== 叠加、分流与推广

- *叠加：* 独立泊松过程 $N_t^(i)$ 的和仍为泊松过程，强度为 $sum_i lambda_i$。
- *分流：* 每次到达独立标记为第 $i$ 类，概率为 $r_i$，则各类计数是相互独立、强度 $lambda r_i$ 的泊松过程。
- *非时齐泊松过程：* 将常强度改为 $lambda(t)>=0$。令 $Lambda(t)=integral_0^t lambda(u) dif u$，则
  $N_t-N_s tilde "Po"(Lambda(t)-Lambda(s))$。
- *复合泊松过程：* 若跳幅 $Y_k$ 独立同分布且与 $N$ 独立，
  $X_t=sum_(k=1)^(N_t)Y_k$。有 $E(X_t)=lambda t E(Y_1)$，
  $D(X_t)=lambda t E(Y_1^2)$。
- *更新过程：* 对一般非负独立同分布间隔 $W_k$，令 $T_n=sum_(k=1)^n W_k$，
  $N_t=max\{n:T_n<=t\}$。指数间隔给出的特殊情形正是泊松过程。
