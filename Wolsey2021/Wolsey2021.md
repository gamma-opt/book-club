Gamma-opt bookclub 2021
===


[![hackmd-github-sync-badge](https://hackmd.io/ZuG05Tv6SSyeKJ6Rdk6j7g/badge)](https://hackmd.io/ZuG05Tv6SSyeKJ6Rdk6j7g)


Notes form the 2021 Gamma-opt book club when we read  Integer programming - Wolsey (2021)

### Participants
- Fabricio Oliveira
- Paula Weller
- Jaan Tollander de Balsch
- Olli Herrala
- Helmi Hankimaa 
- Nikita Belyak


Chapter 1: Formulations
---

## 1.1 Introduction

Classical problems in integer (IP) and mixed-integer (MIP) programming include

- **Scheduling**: sequence activities in time, sometimes allocating resources to it. Can be either with the hope of finding *a* feasible schedule, as in a feasibility problem (train schedule example) or finding the *best* schedule according to a given metric (airline crew example)

- **Production planning**: decide how inputs are converted into outputs, taking into account production system characteristics, such as capacities, minimum batches (production planning and cutting problem examples) and time-related relationships (commitment issues, reserves, and ramping times, as in the electricity generation planning example). 

- **Network design**: design an efficient distribution network in terms of delivery capacity. Might consider aspects related to distribution capacity, node capacity and perhaps combined with some of the above (e.g., supply chain management problems)

- **Assignment, covering and graph related**: problems that can be structured as a graph in which one wishes to identify structures (assignments between nodes, such as in the kidney exchange program example), set covers (radiation therapy and other coverring problems)

Basically, any problem in which decisions are made to be discreted, most notably binary, to represent Yes/ No, True/ False assignments.


## 1.2 What is an integer program?

We depart from the linear program

\begin{equation}
    \text{max. } \left\{cx : Ax \leq b, x \geq 0\right\} \tag{LP}
\end{equation}

where $A$ is a $m \times n$ matrix, $c$ is a $n$-dimensional row vector and $b$ is $m$-dimensional column vector; $x$ is a column vector of (continuous) decision variables

A mixed-integer linear program has the form

\begin{equation}
    \text{max. } \left\{c^\top x + h^\top y : Ax + Gy \leq b, x \geq 0 \text{ and integer}, y \geq 0 \right\} \tag{MIP}
\end{equation}

where $G$ is a $m \times p$ matrix and $h$ is a $p$-dimensional row vector. $x$ is a column vector of $n$ integer decision variables and column vector with $p$ continuous decision variables. 

If $n = 0$, we have a integer programming problem (IP); if furthermore $x \in \{0,1\}$, we have a binary integer programming problem (though everyone would call it an IP anyways).

A **combinatorial optimisation problem** is of the form

\begin{equation}
 \text{min. }_{S \subseteq N} \left\{ \sum_{j \in S} c_j : S \in \mathcal{F} \right\} \tag{COP}
\end{equation}

In COP, the objective is finding a subset from a set of feasible subsets $\mathcal{F}$ from $N$ (i.e, combinations of elements in $N$) such that an attribute of the items $j \in N$ is optimised. 

LP and MIP are very much look-alikes, and that is why LP theory underpins the understanding and solving of MIPs. One idea that immediately comes to place is the notion of rounding, which has its role in the solving of MIPs, but in more sophisticated way than first thought.

![](https://i.imgur.com/mFkHYI0.png) 

![](https://i.imgur.com/4LJmXmr.png)

The optimal integer solution is (5,0), while the optimal LP is that point at the top of the triangle. Notice that no rounding can lead you to the optimal IP solution.


## 1.3 Formulating IPs and BIPs

One powerful tool from MIP is the notion of a $n$-dimensional 0-1 *incidence vector* $x^S$, which is connected to the notion of forming subsets $S \subseteq N$ in COPs:
$$
x_j^S = 1 \Leftrightarrow j \in S; x_j^S = 0 \text{ otherwise.}
$$

Let's look at 4 classic problems in IP:


### The assigment problem

$x_{ij}$ - 1 if person $i$ does job $j$, $x_{ij} = 0$ otherwise.
$c_{ij}$ - assignment cost.

\begin{equation}
    \begin{aligned}
        \text{min. }  & \sum_{i=1}^n \sum_{j=1}^n c_{ij} x_{ij} \\
        \text{s.t.: } & \sum_{j=1}^n x_{ij} = 1, \text{ for } i=1,\dots,n \\
                      & \sum_{i=1}^n x_{ij} = 1, \text{ for } j=1,\dots,n \\
                      & x_{ij} \in \{0,1\}, \text{ for } i=1,\dots,n, j=1,\dots,n.
    \end{aligned} \tag{AP}
\end{equation}


### The 0-1 knapsack problem

$x_{j}$ - 1 if item is selected, $x_{j} = 0$ otherwise.
$c_{j}$ - item value
$b$ - weight capacity

\begin{equation}
    \begin{aligned}
        \text{max. }  & \sum_{j=1}^n c_j x_j \\
        \text{s.t.: } & \sum_{j=1}^n a_j x_j \leq b \\
                      & x_{j} \in \{0,1\}, \text{ for } j=1,\dots,n.
    \end{aligned} \tag{KP}
\end{equation}


### The set covering problem

Set covering is strongly connected with contexts realted to services, including emergencies. Let $M = \{1, \dots, m \}$ be a set of regions and $N = \{ 1, \dots, n\}$ be a set of candidate locations. Let $S_j \subseteq M$ be the regions that are served by a server located at $j \in N$. The equivalent COP is

\begin{equation}
    \text{min. }_{T \subseteq N} \left\{ \sum_{j \in T} c_j : \bigcup_{j \in T} S_j = M \right\} \tag{SCP}
\end{equation}

where $c_j$ is a location cost. As it is the case with all COPs, we can formulate SCP as a 0-1 IP. For that we need an *incidence matrix* $A$ (with elements $a_{ij}$), where $a_{ij} = 1$ if $i \in S_j$, and 0 otherwise. 

\begin{equation}
    \begin{aligned}
        \text{min. }  & \sum_{j=1}^n c_j x_j \\
        \text{s.t.: } & \sum_{j=1}^n a_{ij} x_j \geq 1, \text{ for } i=1,\dots,m\\
                      & x_{j} \in \{0,1\}, \text{ for } j=1,\dots,n.
    \end{aligned} \tag{SCP}
\end{equation}

### The travelling salesman problem

Given a collection of cities $N$, visit each city traversing arcs $(i,j), i,j \in N,$ only once, while minimising the total tour cost 
$$
\sum_{i=1}^n \sum_{j=1}^n c_{ij}.
$$

$x_{ij}$ - 1, if the arc (i,j) is traversed (used to form the tour), 0 otherwise. $x$ is not defined for $i = j$.

\begin{equation}
    \begin{aligned}
        \text{min. }  & \sum_{i=1}^n \sum_{j=1}^n c_{ij} x_{ij} \\
        \text{s.t.: } & \sum_{j=1}^n x_{ij} = 1, \text{ for } i=1,\dots,n \\
                      & \sum_{i=1}^n x_{ij} = 1, \text{ for } j=1,\dots,n \\
                      & x_{ij} \in \{0,1\}, \text{ for } i=1,\dots,n, j=1,\dots,n.
    \end{aligned} \tag{TSP}
\end{equation}

Notice that, as is, this is the same as (AP), which might allow for disconnected subtours between nodes $S \subset N$, as shown below.

![](https://i.imgur.com/H5rq62j.png)

Either of the additional constraints can be used to prevent this: 

1. **The cut-set constraint**
$$
\sum_{i \in S} \sum_{j \not\in S} x_{ij} \geq 1 \text{ for } S \subset N, S \not= \emptyset
$$

2. **The subtour elimination constraint**
$$
\sum_{i \in S} \sum_{j \in S} x_{ij} \leq |S|-1 \text{ for } S \subset N, 1 < |S| < n.
$$
 
 
## 1.4 The combinatorial explosion

Feasible space of some of the above problems:
- AP: $n!$
- KP and SCP: $2^n$
- TSP: $(n-1)!$

![](https://i.imgur.com/ux6xbqe.png)

when I think about those, I like to contrast these numbers with these: https://www.physicsoftheuniverse.com/numbers.html


## 1.5 Mixed integer formulations

### Modelling fixed costs

Binary variables can be used to model specific nonlinear "behaviours", not necessarily representing an actual decision. For example, the function $h(y)$ 

![](https://i.imgur.com/8Z7Lhx1.png)

can be modelled as $fx + py$ with the additional constraints $y \leq Cx, x \in \{0,1\}$. In this case specifically, notice that under an optimisation standpoint, it makes no sense to have $y=0$ while $x=1$, though in principle it would be a feasible solution. Otherwise, a more sophisticated set of constraints (still IP) would be needed.


### Uncapacited facility location

MIP combining netork desing and production planning aspects. 

- $x_j = 1$, if depot $j \in N$ is used, 0 otherwise.
- $y_{ij}$ - fraction of the demand of client $i \in M$ satisfied from depot $j$.

\begin{equation}
    \begin{aligned}
        \text{min. }  & \sum_{i \in M} \sum_{j \in N} c_{ij} y_{ij} + \sum_{j \in N}f_j x_j \\
        \text{s.t.: } & \sum_{i \in M} y_{ij} \leq m x_j, \text{ for } i=1,\dots,m \\
                      & \sum_{j=1}^n y_{ij} = 1, \text{ for } i=1,\dots,m \\
                      & x_{j} \in \{0,1\}, \text{ for } j=1,\dots,n \\
                      & y_{ij} \ge 0, \text{ for } i=1,\dots,n, j=1,\dots,n.
    \end{aligned} \tag{UFL}
\end{equation} 


### Uncapacited lot sizing

Lot sizing is the common nomenclature for a classical production planning setting. The input data is
- $f_t$ - fixed cost of producing in time period $t \in \{1, \dots, n\}$
- $p_t$ - unit production cost in time period $t \in \{1, \dots, n\}$
- $h_t$ - unit storage cost in time period $t \in \{1, \dots, n\}$
- $d_t$ - demand in time period $t \in \{1, \dots, n\}$

The decision variables are:
- $y_t$ is the amount produced in period $t$
- $s_t$ is the amount produced in period $t$
- $x_t = 1$, if production occurs in $t$, 0 otherwise.

\begin{equation}
    \begin{aligned}
        \text{min. }  & \sum_{t=1}^n p_t y_t + \sum_{t=1}^n h_t s_t + \sum_{t=1}^n f_t x_t \\
        \text{s.t.: } & s_{t-1} + y_t = d_t + s_t, \text{ for } t = 1, \dots, n \\
                      & y_{t} \leq Mx_t, \text{ for } t = 1, \dots, n \\
                      & s_t, y_t \geq 0, \text{ for } t = 1, \dots, n  \\
                      & x_t \in \{0,1\}, \text{ for } t = 1, \dots, n.
    \end{aligned} \tag{ULS}
\end{equation} 

Notice that in the absence of an upper bound for $y_{t}$ a sufficiently large constant $M$ must be used, which can lead to problems related numerical issues if not carefully considered.

### Discrete alternatives or disjuctions 

Useful for modelling nonconvex affine regions. Assume that either
$$ 
a^1 y\le b_1 \text { or } a^2 y\le b_2 
$$
must hold, with $0 \leq y \leq u$. Such case is illustrated in the picture below.

![](https://i.imgur.com/Pc0v8Ef.png)

Define $M \geq \max\{a^i - b_i : 0 \leq y \leq u\}$ for $i=1,2$. Then, we can use the constraints
\begin{equation}
    \begin{aligned}
        & a^iy - b_i \leq M(1 - x_i), \text{ for } i=1,2 \\
        & x_1 + x_2 = 1 \\
        & x_i \in \{0,1\}, , \text{ for } i=1,2.
    \end{aligned} 
\end{equation}

This type of disjunction arise in scheduling problems, for example. Let $t_i$ be a starting time for processing a job $i$ and $p_i$ the processing time. Then either
$$
t_2 \geq t_1 + p_1 \text{ or } t_1 \geq p_2 + t_2.
$$


## 1.6 Alternative formulations

Typically, an IP has a multitude of possibe correct formulations. Let us look into what makes a formulation better than another. The definitions below state the notion of a polyhedral set (polyhedron) and exaclty what is a *formulation. 

![](https://i.imgur.com/zU285s2.png)

**Example 1.2**: Let $X = \{(1,1), (2,1), (3,1), (1,2), (2,2), (3,2), (2,3) \}$. The picture show alternative formulations for $X$.

![](https://i.imgur.com/3fVPtA0.png)


### Equivalent formulation for UFL

This is an example in which the alternative formulations use the same variables. 

Notice that the constraint
$$
\sum_{i \in M} y_{ij} \leq m x_j,  \text{ for } j \in N 
$$
can be equivalently expressed as
$$
0 \le y_{ij} \le x_j, \text{ for } i \in M, j \in N 
$$
leading to the alternative formulation
\begin{equation}
    \begin{aligned}
        \text{min. }  & \sum_{i \in M} \sum_{j \in N} c_{ij} y_{ij} + \sum_{j \in N}f_j x_j \\
        \text{s.t.: } & y_{ij} \leq  x_j, \text{ for } i \in M, j \in N\\
                      & \sum_{j=1}^n y_{ij} = 1, \text{ for } i \in M \\
                      & x_{j} \in \{0,1\}, \text{ for } j \in N \\
                      & y_{ij} \ge 0, \text{ for } i \in M, j \in N.
    \end{aligned} \tag{UFL}
\end{equation} 


### Extended formulation for ULS

Theer are cases in which alternative formulations are obtained by means of different decision variables, which makes the comparison between them somewhat tricky. These are called extended formulations. 

Let us redefine the variables in the ULS:
- $w_{it}$ - amount produced in period $i$ to satisfy the demand at period $t$.
- $x_t$ - if production occurs in period $t$ (like before).

And redefine the parameter $c_t = p_t + h_t + \dots + h_n$. Then the ULS problem can be reformulated as

\begin{equation}
    \begin{aligned}
        \text{min. }  & \sum_{i=1}^n \sum_{t=1}^n c_i w_{it} + \sum_{t=1}^n f_t x_t \\
        \text{s.t.: } & \sum_{i=1}^t w_{it} = d_t, \text{ for } t = 1, \dots, n \\
                      & w_{it} \leq d_tx_t, \text{ for } i \leq j \text { and }i,t = 1, \dots, n \\
                      & w_{ij} \geq 0, \text{ for } i \leq j \text { and }i,t = 1, \dots, n \\
                      & x_t \in \{0,1\}, \text{ for } t = 1, \dots, n.
    \end{aligned} \tag{ULS}
\end{equation} 


## 1.7 Good and ideal formulation

The answer for understanding how a formulation can be better than another comes from geometrical insight.

Ideally, we want to be able to solve a MIP problem by solving its LP formulation, $P$, since solving LPs is "easy" (certainly easier than MIPs) from a complexity and computational standpoint. The picture below show what a ideal formulation $P$ would look like:

![](https://i.imgur.com/eiBzcgK.png)

This formalises the notion of ideal formulations:

![](https://i.imgur.com/qbfIeSz.png)

**Proposition 1.1** can be seen by noticing that, by **Definition 1.3**, convex hulls are polyhedral sets in general. **Proposition 1.2** is a central result from linear optimisation. 

Combining **Proposition 1.1** and **Proposition 1.2** leads us to the conclusion that we can replace the (M)IP
$$
\text{max.}\{cx : x \in X\} \text{ with } \text{max.}\{cx : x \in conv(X)\},
$$
which would allow us to solve the MIP as an LP. However, unless we already have an ideal formulation at hand (meaning that $P = conv(X)$),  this is only of theoretical value, as describing $conv(X)$ require a (typically too) large number of constraints that are not trivially identifiable.

But not everything is lost. We can use the above idea to define how one formulation can be better than another.

![](https://i.imgur.com/YwGCynb.png)

### Formulations for the UFL

Let $P_1$ be the formulation for the UFL with the aggregate constraint
\begin{equation}
    \sum_{i \in M} y_{ij} \leq m x_j, \tag{1}
\end{equation}
and $P_2$ be the formulation with the $m$ individual constraints
\begin{equation}
    y_{ij} \leq x_j, \text{ for }i \in M. \tag{2}
\end{equation}

The reasoning of $P_2 \subset P_1$ (i.e., that $P_2$ is better than $P_1$) is as follows.
1. Any point $(x,y)$ that satisfies $y_{ij} \leq x_j, \text{ for }i \in M$, also satisfies $\sum_{i \in M} y_{ij} \leq m x_j$ by simply summing in $i \in M$. Thus $P_2 \subseteq P_1$
2. Now we need to show that there are points in $P_1$ that are not in $P_2$. For that, let (x,y) be such that:
    - Suppose that $n$ divides $m$: $m = k \times n$ with $k > 1$.
    - For example: $m=4, n=2$. Assemble the solution: $x_j = k/m$ and $y_{ij}$ as
        \begin{equation}
            y_{11}, y_{21}, y_{32}, y_{42} = 1, 0 \text{ otherwise.}
        \end{equation}
    - This solution satisfy (1) but does not satisfy (2). Thus, $(x,y) \in P_1 \setminus P_2$, implying that $P_2 \subset P_1$.     
        

### Formulations for the ULS

In the case of the ULS, things are more complicated because the problems have different variables. The way to go around that is to rely on *projections*.

![](https://i.imgur.com/G4cAvo5.png)

Let $P_1$ be defined as 
\begin{equation}
    \begin{aligned}
        & s_{t-1} + y_t = d_t + s_t, \text{ for } t = 1, \dots, n \\
        & y_{t} \leq Mx_t, \text{ for } t = 1, \dots, n \\
        & s_t, y_t \geq 0, \text{ for } t = 1, \dots, n \\
        & 0 \leq x_t \leq 1, \text{ for } t = 1, \dots, n.
    \end{aligned} \tag{$P_1$}
\end{equation} 

Recall that $M = \sum_{t} d_t$, which is important in the argument made later on. Now $P_2$ would need to be the projection of $Q_2$ onto variables $(y,s,x)$ so we can compare them. $Q_2$ can be defined as:

\begin{equation}
    \begin{aligned}
        & s_{t-1} + y_t = d_t + s_t, \text{ for } t = 1, \dots, n \\
        & \sum_{i=1}^t w_{ij} = d_t, \text{ for } t = 1, \dots, n \\
        & w_{it} \leq d_i x_i, , \text{ for } i \leq j \text { and }i,t = 1, \dots, n \\
        & y_i = \sum_{t=i}^n w_{it}, \text{ for } t = 1, \dots, n \\
        & w_{it} \ge 0 \text{ for } i \leq j \text { and }i,t = 1, \dots, n \\
        & 0 \leq x_t \leq 1,, s_t \ge 0, y_y \ge 0, \text{ for } t = 1, \dots, n.
    \end{aligned} \tag{$Q_2$}
\end{equation} 

One way to see the projection of ($Q_2$) onto $(y,s,x)$ is to simply think of $w_{it}$ as "fixed".

Considering the point $(\overline{x},\overline{s}, \overline{y})$ with $\overline{y}_t = d_t$, $\overline{x}_t = d_t/M$, we can see that it satisfies all constraints in $P_1$. 

Now, if you sum in $w_{it} \leq d_i x_i$ in $t$, and substitute $(\overline{x},\overline{s}, \overline{y})$ you obtain

\begin{equation}
    \sum_{i=1}^n w_{it} \leq \sum_{i=1}^n d_i \overline{x}_i,\text{ for } t = 1, \dots, n.
\end{equation}

By taking $t < n$, we have 
\begin{equation}
    \sum_{t=1}^n w_{it} \leq d_t \frac{\sum_{i=1}^n d_i}{M} < d_t
\end{equation}

which violates the second constraint $\sum_{i=1}^t w_{ij} = d_t$ in $P_2$. To see that a point that satisfies $P_2$ also satisfy $P_1$, you can sum the third constraint in $Q_2$ in $t$, obtaining 

$$
\sum_{i=1}^n w_{it} = \overline{y}_i \leq \sum_{i=1}^n d_i \overline{x}_i = M\overline{x}_i,\text{ for } t = 1, \dots, n.
$$
