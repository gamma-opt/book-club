Gamma-opt bookclub 2021
===


Notes form the 2021 Gamma-opt book club when we read  Integer programming - Wolsey (2021)

### Participants
- Fabricio Oliveira
- Paula Weller
- Jaan Tollander de Balsch
- Olli Heralla
- Helmi Hankimaa
- Nikita Belyak


Chapter 1: Formulations
---

## 1.1 Introduction

Classical problems in integer (IP) and mixed-integer (MIP) programming include

- **Scheduling**: sequence activities in time, sometimes allocating resources to it. Can be either with the hope of finding *a* feasible schedule, as in a feasibility problem (train schedule example) or finding the *best* schedule according to a given metric (airline crew example)

- **Production planning**: decide how inputs are converted into outputs, taking into account production system characteristics, such as capacities, minimum batches (production planning and cutting problem examples) and time-related relationships (commitment issues, reserves, and ramping times, as in the electricity generation planning example). 

- **Network design**: design an efficient distribution network in terms of delivery capacity. Might consider aspects related to distribution capacity, node capacity and perhas combined with some of the above (e.g., suplly chain management problems)

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
    \text{max. } \left\{c^\top x + h^\top y : Ax + Gy \leq b, x \geq 0 \text{and integer}, y \geq 0 \right\} \tag{MIP}
\end{equation}

where $G$ is a $m \times p$ matrix and $h$ is a $p$-dimensional row vector. $x$ is a column vector of $n$ integer decision variables and column vector with $p$ continuous decision variables. 

If $n = 0$, we have a integer programming problem (IP); if furthermore $x \in \{0,1\}$, we have a binary integer programming problem (though everyone would call it an IP anyways).

A **combinatorial optimisation problem** is of the form

\begin{equation}
 \text{min. }_{S \subseteq N} \left\{ \sum_{j \in S} c_j : S \in \mathcal{F} \right\} \tag{COP}
\end{equation}

In COP, the objective is finding a subset from a set of feasible subsets $\mathcal{F}$ from $N$ (i.e, combinations of elements in $N$) such that an attribute of the items $j \in N$ is optimised. 

LP and MIP are very much look-alikes, and that is whyn LP theory underpins the understanding and solving of MIPs. One idea that immediately comes to place is the notion of rounding, which has its role in the solving of MIPs, but in more sophisticated way than first thought.

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

Set coverin is strongly connected with contexts realted to services, including emergencies. Let $M = \{1, \dots, m \}$ be a set of regions and $N = \{ 1, \dots, n\}$ be a set of candidate locations. Let $S_j \subseteq M$ be the regions that are served by a server located at $j \in N$. The equivalent COP is

\begin{equation}
    \text{min. }_{T \subseteq N} \left\{ \sum_{j \in T} c_j : \bigcup_{j \in T} S_j = M \right\} \tag{SCP}
\end{equation}

where $c_j$ is a location cost. As it is the case with all COPs, we can formulate SCP as a 0-1 IP. For that we need an *incidence matrix* $A$ (with elements $a_{ij}$), where $a_{ij} = 1$ if $i \in S_j$, and 0 otherwise. 

\begin{equation}
    \begin{aligned}
        \text{min. }  & \sum_{j=1}^n c_j x_j \\
        \text{s.t.: } & \sum_{j=1}^n a_{ij} x_j \geq 1, \text{ for } i=1,\dots,n\\
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

can be modelled as $fx + py$ witht he additional constraints $y \leq Cx, x \in \{0,1\}$. In this case specifically, notice that under an optimisation standpoint, it makes no sense to have $y=0$ while $x=1$, though in principle it would be a feasible solution. Otherwise, a more sophisticated set of constraints (still IP) would be needed.


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
    \end{aligned} \tag{UFL}
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
\sum_{i \in M} y_{ij},  \text{ for } j \in N 
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

