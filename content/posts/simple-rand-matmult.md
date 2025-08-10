---
title: "A Simple Randomized Matrix Multiplication Algorithm 1"
date: "2025-08-09"
tags: []
categories: []
ShowToc: false
TocOpen: false
math: true
---

> 千里之行，始於足下。
>
> [The journey of a thousand miles begins with a simple step.]
>
> -- *Lao Tzu*

# Introduction

We all have to start somewhere, and I have decided to start the blog with something simple, that most people in the quantitative sciences would have to deal with at some point in their lives: matrix multiplication.

Performing matrix multiplication fast is important for many applications, like machine learning or numerical methods. We will consider matrices (square of order $n$, for simplicity) $A=(a_{ij})\_{i,j=1}^n$, and $B=(b_{ij})\_{i,j=1}^n$, think about ways to evaluate $AB$. The simplest, most naive algorithm would be as follows:

***
1. Initialize an $n\times n$ zero matrix $C=(c_{ij})_{i,j=1}^n$.
2. * For $i=1,\dots,n$,
        * For $j=1,\dots,n$,
            * For $k=1,...,n$,
                * $c_{ij} \mathrel{+}= a_{ik}b_{kj}$.
3. Return $C$.
***

Naturally, because of the triple for-loop in the algorithm, this runs in $\mathcal{O}(n^3)$ time. This can be considered abysmally slow! For many years, people could not come up with faster ways to multiply matrices, but Volker Strassen published [his algorithm](https://en.wikipedia.org/wiki/Strassen_algorithm) with a runtime of $\mathcal{O}(n^{\log_2 7})$ in 1969, and demonstrated that the $\mathcal{O}(n^3)$ barrier could indeed be broken. Nowadays, several more faster algorithms have been discovered, with the [best bound](https://arxiv.org/abs/2404.16349) (as of writing) being $\mathcal{O}(n^{2.371339})$.

But what if this is not enough? What if we wanted *faster*? 

Generally speaking, one can judge an algorithm in three ways: speed, simplicity, and correctness. Unfortunately, oftentimes these three goals stand in contradiction to each other. Here, we need more speed, and people have already racked their brains hard and came up with the most complicated fast algorithms for matrix multiplication. What we will do today is sacrifice a little *correctness*. More precisely, we will introduce some randomization into the picture, and only ask that our algorithm returns something (quantifiably) close to the right answer most of the time.

# The Idea
Of course, we all know how to multiply matrices, so observe that

$$
    AB = \left(\sum_{k=1}^n a_{ik}b_{kj}\right)\_{i,j=1}^n = \sum_{k=1}^n\left(a\_{ik}b_{kj}\right)_{i,j=1}^n.
$$
The right-hand side is a product of rank-one matrices. The idea is that we could randomly sample the rank-one matrices (essentially, we are randomly sampling $k$ from $1,\dots,n$) and sum them up. Evaluating each rank-one matrix takes $\mathcal{O}(n^2)$ time, so if we can show that we do not have to sample too many times, we will have a fast algorithm. Therefore, we will consider the following algorithm template.


***
Inputs: the sample size $S\in\mathbb{Z}_{>0}$, and probabilities $p_1,...,p_n>0$ summing to $1$.
1. Initialize an $n\times n$ zero matrix $C=(c_{ij})_{i,j=1}^n$.
2. * Repeat $S$ times:
        * Pick a random $k\in\{1,\dots,n\}$ with probability $p_k$.
        * For $i=1,\dots,n$,
            * For $j=1,...,n$,
                * $c_{ij} \mathrel{+}= \frac{1}{p_i} a_{ik}b_{kj}$.
3. Return $D := \frac{1}{S}C$.
***

Naturally, one has to pick the $p_k$ so that the algorithm works. But by magic, the choice of $p_k$ does not matter in the following sense:

**Lemma.** Regardless of the choice of $p_k$, one has $\mathbb{E}[D] = AB$.

In English, this means that $D$ is always an unbiased estimator of $AB$.

*Proof.* We let $k_m$ be the $m$'th sample of $k$, so that

$$
\begin{align*}
D=\frac{1}{S}\sum_{m=1}^S\frac{1}{p_{k_m}}(a_{ik_m}b_{k_mj})_{i,j=1}^n.
\end{align*}
$$

By linearity of expectation, one then gets

$$
\begin{align*}
\mathbb{E}[D]
&=\frac{1}{S}\sum_{m=1}^S\mathbb{E}\left[\frac{1}{p_{k_m}}(a_{ik_m}b_{k_mj})\_{i,j=1}^n\right] \\\\
&=\frac{1}{S}\sum_{m=1}^S\sum_{k=1}^n p_k \frac{1}{p_{k}}(a_{ik}b_{kj})\_{i,j=1}^n \\\\
&=\frac{1}{S}\sum_{m=1}^S\sum_{k=1}^n (a_{ik}b_{kj})\_{i,j=1}^n \\\\
&=\frac{1}{S}\sum_{m=1}^S AB \\\\
&= AB.
\end{align*}
$$

So our algorithm technically works. Unfortunately, works is not enough: the variance of the outputs could be absurdly high, such that most of the answers are way off from the correct answer. The next thing to do is to control this variance. This means we got to control $D-AB$, which also means we have to impose a notion of size onto $D-AB$. Our choice is to use the Frobenius norm:

$$
\\|A\\|\_{F} = \sqrt{\sum_{i,j=1}^n |a_{ij}|^2}.
$$

One has, by the independence of the $k_m$, that
$$
\begin{align*}
\mathbb{E}[\\|D-AB\\|\_F^2] 
&= \sum_{i,j=1}^n \mathbb{E}\left[\left(\frac{1}{S}\sum_{m=1}^S \frac{1}{p_{k_m}}a_{ik_m}b_{k_mj} - \sum_{k=1}^na_{ik}b_{kj}\right)^2\right] \\\\
&= \sum_{i,j=1}^n \mathrm{Var}\left[\frac{1}{S}\sum_{m=1}^S \frac{1}{p_{k_m}}a_{ik_m}b_{k_mj}\right] \\\\
&= \sum_{i,j=1}^n \frac{1}{S^2}\sum_{m=1}^S\mathrm{Var}\left[\frac{1}{p_{k_m}}a_{ik_m}b_{k_mj}\right] \\\\
&= \sum_{i,j=1}^n \frac{1}{S^2}\sum_{m=1}^S\sum_{k=1}^n\left(\frac{a_{ik}^2b_{kj}^2}{p_k} - a_{ik}^2b_{kj}^2\right) \\\\
&= \sum_{i,j=1}^n \frac{1}{S}\sum_{k=1}^n\left(\frac{a_{ik}^2b_{kj}^2}{p_k} - a_{ik}^2b_{kj}^2\right) \\\\
&= \frac{1}{S}\sum_{k=1}^n\frac{1}{p_k}\left(\sum_{i=1}^n a_{ik}^2\right)\left(\sum_{j=1}^n b_{kj}^2\right) - \frac{1}{S}\\|AB\\|_F.
\end{align*}
$$

It is in our interest to minimize this expression. Fortunately, one can easily lower bound the expression via a clever use of Cauchy-Schwarz:

$$
\begin{align*}
&\frac{1}{S}\sum_{k=1}^n\frac{1}{p_k}\left(\sum_{i=1}^n a_{ik}^2\right)\left(\sum_{j=1}^n b_{kj}^2\right) - \frac{1}{S}\\|AB\\|\_F \\\\
&= \frac{1}{S}\left(\sum_{k=1}^n p_k\right)\left(\sum_{k=1}^n\frac{1}{p_k}\left(\sum_{i=1}^n a_{ik}^2\right)\left(\sum_{j=1}^n b_{kj}^2\right)\right) - \frac{1}{S}\\|AB\\|\_F \\\\
&\geq \frac{1}{S}\sum_{k=1}^n\sqrt{\sum_{i=1}^n a_{ik}^2}\sqrt{\sum_{j=1}^n b_{kj}^2} - \frac{1}{S}\\|AB\\|\_F,
\end{align*}
$$
then recall that the equality in Cauchy-Schwarz holds if and only if we set $p_k$ proportional to $\sqrt{\sum_{i=1}^n a_{ik}^2}\sqrt{\sum_{j=1}^n b_{kj}^2}$. Therefore, we minimize the expected error if we set
$$
p_k = \frac{\sqrt{\sum_{i=1}^n a_{ik}^2}\sqrt{\sum_{j=1}^n b_{kj}^2}}{\sum_{k'=1}^n\sqrt{\sum_{i=1}^n a_{ik'}^2}\sqrt{\sum_{j=1}^n b_{k'j}^2}}.
$$

The $p_k$ are known as the *optimal sampling probabilities*.

Algorithms like ours tend to be rather resilient to some variation in the sampling probabilities. To illustrate this point in the next post ([here](/personalpage/posts/simple-rand-matmult-2/)), we will consider instead *$\beta$-approximately optimal sampling probabilities* for some $0<\beta<1$. That is, we will be interested in sampling probabilities $p_k$ satisfying
$$
p_k \geq \frac{\beta\sqrt{\sum_{i=1}^n a_{ik}^2}\sqrt{\sum_{j=1}^n b_{kj}^2}}{\sum_{k'=1}^n\sqrt{\sum_{i=1}^n a_{ik'}^2}\sqrt{\sum_{j=1}^n b_{k'j}^2}}.
$$
Under this assumption, we can quickly derive, using Cauchy-Schwarz,
$$
\begin{align*}
\mathbb{E}[\\|D-AB\\|\_F^2] 
&= \frac{1}{S}\sum_{k=1}^n\frac{1}{p_k}\left(\sum_{i=1}^n a_{ik}^2\right)\left(\sum_{j=1}^n b_{kj}^2\right) - \frac{1}{S}\\|AB\\|\_F \\\\
&\leq \frac{1}{S}\sum_{k=1}^n\frac{1}{p_k}\left(\sum_{i=1}^n a_{ik}^2\right)\left(\sum_{j=1}^n b_{kj}^2\right) \\\\
&\leq \frac{1}{\beta S}\sum_{k=1}^n\sqrt{\sum_{i=1}^n a_{ik}^2}\sqrt{\sum_{j=1}^n b_{kj}^2}\sum_{k'=1}^n \sqrt{\sum_{i=1}^n a_{ik'}^2}\sqrt{\sum_{j=1}^n b_{k'j}^2} \\\\
&= \frac{1}{\beta S}\left(\sum_{k=1}^n\sqrt{\sum_{i=1}^n a_{ik}^2}\sqrt{\sum_{j=1}^n b_{kj}^2}\right)^2 \\\\
&\leq \frac{1}{\beta S}\left(\sum_{k=1}^n\sum_{i=1}^n a_{ik}^2\right)\left(\sum_{k=1}^n\sum_{j=1}^n b_{kj}^2\right) \\\\
&= \frac{1}{\beta S}\\|A\\|_F^2\\|B\\|_F^2.
\end{align*}
$$

# What Next?

So far, choosing the $p_k$ optimally has helped us to minimize the expected error. But we have not managed to control the variance of the error yet. In other words, we want the error to be low, with high probability. This is where picking the value of $S$ comes in. The act of collecting enough samples to kill the variance is usually called *concentration* in the literature, and this is typically established by invoking a concentration inequality such as Markov's inequality, Chebyshev's inequality, or some Chernoff/Hoeffding bound. In the next point, we will introduce McDiarmid's inequality, and using it we will perform concentration to get a good value of $S$.

# References
Micheal Mahoney's Lecture Notes on Randomized Linear Algebra, over at <https://arxiv.org/pdf/1608.04481>.