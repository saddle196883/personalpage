---
title: "A Simple Randomized Matrix Multiplication Algorithm 2"
date: "2025-08-10"
tags: []
categories: []
ShowToc: false
TocOpen: false
math: true
---

> The huger the mob, and the greater the apparent anarchy, the more perfect is its sway. It is the supreme law of Unreason. 
>
> -- *Sir Francis Galton*

This post is the continuation of the previous post over [here](/personalpage/posts/simple-rand-matmult/).

# Last Time
In the previous post, we began the analysis of a randomized matrix multiplication algorithm. We managed to prove that the output of the algorithm gives us our product $AB$ in expectation, and furthermore, that the expected squared deviation $\mathbb{E}[\\|D-AB\\|_F^2]$ satisfies

$$
\begin{align*}
\mathbb{E}[\\|D-AB\\|\_F^2] \leq \frac{1}{\beta S}\\|A\\|_F^2\\|B\\|_F^2,
\end{align*}
$$
when we sample with $\beta$-optimal sampling probabilities, i.e.

$$
p_k \geq \frac{\beta\sqrt{\sum_{i=1}^n a_{ik}^2}\sqrt{\sum_{j=1}^n b_{kj}^2}}{\sum_{k'=1}^n\sqrt{\sum_{i=1}^n a_{ik'}^2}\sqrt{\sum_{j=1}^n b_{k'j}^2}}.
$$

Now that we know $\mathbb{E}[\\|D-AB\\|\_F^2]$ in expectation, our next step is to ensure that $\\|D-AB\\|\_F^2$ does not deviate too far from its expected value most of the time. In order to do so, we will need some concentration inequality.

# McDiarmid's inequality

McDiarmid's inequality is typically used to control the deviations of functions from their expected values, if the function does not "vary too much" if we vary one variable while leaving the others fixed. More concretely, we need our function to satisfy the following property:

**Definition.** A function $f:\mathcal{X}\_1\times \dots \mathcal{X}\_n\to \mathbb{R}$ satisfies the *bounded difference property* if and only if there exists constants $c_1,\dots, c_n$ such that for every $x_1\in \mathcal{X}\_1,\dots,x_n\in \mathcal{X}\_n$ and any $i=1,\dots,n$, one has
$$
\sup_{y\in X_i}|f(x_1,\dots,x_i,\dots,x_n) - f(x_1,\dots,y,\dots,x_n)|\leq c_i.
$$

Then McDiarmid's inequality is as follows:

**Theorem. (McDiarmid's inequality)** Let $X_1,\dots,X_n$ be independent random variables taking values in $\mathcal{X}_1,\dots,\mathcal{X}_n$ respectively, and let $f:\mathcal{X}_1\times\dots\times \mathcal{X}_n\to\mathbb{R}$ satisfy the bounded difference property. Then

$$
\mathrm{Pr}\left[f(X_1,\dots,X_n)-\mathbb{E}[f(X_1,\dots,X_n)]\geq \tau\right]\leq \exp\left(-\frac{2\tau^2}{\sum_{i=1}^n c_i^2}\right).
$$
We omit the proof, but the typical way to prove McDiarmid's inequality is to construct a Doob martingale and apply the Azuma-Hoeffding inequality. (In fact, the original plan for this post was to do precisely that, but I gave up writing about martingales.) 

Concentration inequalities like this one control the probability that one ends up in the tail ends of a distribution. There is a plethora of concentration inequalities out there: check out the [wikipedia page](https://en.wikipedia.org/wiki/Concentration_inequality) for a few of them.

# Concentration

Having stated the inequality, we should now consider how to use it. The value we want to control is $\\|D-AB\\|\_F$. This value depends on our choice of $k_m$ for each $m=1,\dots,S$. It therefore makes sense to define $f:\\{1,...,n\\}^S\to \mathbb{R}$ by
$$
f(k_1, \dots, k_S) = \\|D - AB\\|\_F,
$$
where $D=\sum_{m=1}^S \frac{1}{p_{k_m}}(a_{ik_m}b_{k_mj})_{i,j=1}^n$.

Now suppose that for some fixed $m$, we replaced $k_m$ with $k_m\'$, and called this modified matrix $D\'$. Then by the triangle inequality and Cauchy-Schwarz, we get:
$$
\begin{align*}
\\|D-D\'\\|\_F
&= \left\\|\frac{1}{Sp_{k_m}}(a_{ik_m}b_{k_mj})\_{i,j=1}^n - \frac{1}{Sp_{k_m\'}}(a_{ik_m\'}b_{k_m\'j})\_{i,j=1}^n\right\\|\_F \\\\
&\leq \frac{1}{Sp_{k_m}}\\|(a_{ik_m}b_{k_mj})\_{i,j=1}^n\\| + \frac{1}{Sp_{k_m\'}}\\|(a_{ik_m\'}b_{k_m\'j})\_{i,j=1}^n\\|\_F \\\\
&= \frac{1}{Sp_{k_m}}\sqrt{\sum_{i,j=1}^n a_{ik_m}^2b_{k_mj}^2} + \frac{1}{Sp_{k_m\'}}\sqrt{\sum_{i,j=1}^na_{ik_m\'}^2b_{k_m\'j}^2} \\\\
&= \frac{2}{S} \max_{k=1,\dots,n} \frac{\sqrt{\sum_{i,j=1}^n a_{ik}^2b_{kj}^2}}{p_k} \\\\
&= \frac{2}{S} \max_{k=1,\dots,n} \frac{\sqrt{\sum_{i=1}^n a_{ik}^2}\sqrt{\sum_{j=1}^nb_{kj}^2}}{p_k} \\\\
&\leq \frac{2}{\beta S} \sum_{k=1}^n \sqrt{\sum_{i=1}^n a_{ik}^2}\sqrt{\sum_{j=1}^nb_{kj}^2} \\\\
&\leq \frac{2}{\beta S} \sqrt{\sum_{i,k=1}^n a_{ik}^2}\sqrt{\sum_{j,k=1}^nb_{kj}^2} \\\\
&\leq \frac{2}{\beta S} \\|A\\|\_F\\|B\\|\_F.
\end{align*}
$$
Consequently,
$$
\\|D-AB\\|\_F \leq \\|D\'-AB\\|\_F + \\|D-D\'\\|\_F,
$$
$$
\\|D\'-AB\\|\_F \leq \\|D-AB\\|\_F + \\|D-D\'\\|\_F,
$$
give
$$
\max_{y=1,\dots,n}|f(k_1,\dots,k_m,\dots,k_S) - f(k_1,\dots,y,\dots,k_S)|\leq \frac{2}{\beta S}\\|A\\|_F\\|B\\|_F.
$$

Hence by McDiarmid's inequality one has, for all $\varepsilon>0$,
$$
\mathrm{Pr}\left[\\|D-AB\\|\_F-\mathbb{E}[\\|D-AB\\|\_F]\geq \varepsilon\\|A\\|\_F\\|B\\|\_F\right]\leq \exp\left(-\frac{\varepsilon^2\beta^2S}{2}\right).
$$

On the other hand, because $\mathbb{E}[\\|D-AB\\|_F]\leq \sqrt{\mathbb{E}[\\|D-AB\\|_F^2} \leq\frac{1}{\sqrt{\beta S}}\\|A\\|\_F\\|B\\|\_F$ by Jensen's inequality,
$$
\begin{align*}
&\mathrm{Pr}\left[\\|D-AB\\|\_F-\mathbb{E}[\\|D-AB\\|\_F]\geq \varepsilon\\|A\\|\_F\\|B\\|\_F\right] \\\\
&\geq \mathrm{Pr}\left[\\|D-AB\\|\_F - \frac{1}{\beta S}\\|A\\|\_F\\|B\\|\_F\geq \varepsilon\\|A\\|\_F\\|B\\|\_F\right] \\\\
&= \mathrm{Pr}\left[\\|D-AB\\|\_F \geq \left(\frac{1}{\sqrt{\beta S}} + \varepsilon\right)\\|A\\|\_F\\|B\\|\_F\right].
\end{align*}
$$

Therefore, for $0<\delta<1$, if we set
$$S=\max\left\\{\varepsilon^{-2}\beta^{-1}, 2\varepsilon^{-2}\beta^{-2}\log\frac{1}{\delta}\right\\} = 2\varepsilon^{-2}\beta^{-2}\log\frac{1}{\delta},$$
we will attain
$$
\begin{align*}
\mathrm{Pr}\left[\\|D-AB\\|\_F \geq 2\varepsilon\\|A\\|\_F\\|B\\|\_F\right] \leq \delta.
\end{align*}
$$

In particular, we observe that if the value of $\beta$ does not deviate far from $1$, then its impact on the performance of the algorithm is not too large.
# Summary

Our final algorithm is therefore as follows:
***
1. Compute $w_k=\sqrt{\sum_{i=1}^n a_{ik}^2}\sqrt{\sum_{j=1}^n b_{kj}^2}$ for each $k=1,\dots,n$.
2. Compute $p_k = \frac{w_k}{\sum_{k\'=1}^n w_k\'}$.
2. Initialize an $n\times n$ zero matrix $C=(c_{ij})_{i,j=1}^n$.
3. * Repeat $S=2\varepsilon^{-2}\log\frac{1}{\delta}$ times:
        * Pick a random $k\in\{1,\dots,n\}$ with probability $p_k$.
        * For $i=1,\dots,n$,
            * For $j=1,...,n$,
                * $c_{ij} \mathrel{+}= \frac{1}{p_i} a_{ik}b_{kj}$.
4. Return $D := \frac{1}{S}C$.
***

And from above, we have the following result:

**Theorem.** With probability at least $1-\delta$, the output $D$ of the algorithm above satisfies $\\|D-AB\\|_F < 2\varepsilon\\|A\\|\_F\\|B\\|_F$.

# Final Remarks
Because each $w_k$ is computed in  $\mathcal{O}(n)$, the final algorithm has runtime $\mathcal{O}(\varepsilon^{-2}n^2\log\frac{1}{\delta})$. While on paper this might look nice, it is not very useful in practice: the  $\varepsilon^{-2}$ dooms the algorithm to produce low accuracy results, unless we increase the runtime significantly. This $\varepsilon^{-2}$ is present due to the central limit theorem, and in general cannot be avoided: such is the limitation of pretty much all sampling algorithms out there. Nonetheless, it is a good window into the myriad of possible cuts in runtime that randomization can being into matrix algorithms.

I should remark at the end that at some point in the computations, I had deviated signifcantly from the references that I have been using. In particular, if there are any errors in the computations, it is solely my fault. Please contact me if you do spot any errors.

# References
Micheal Mahoney's Lecture Notes on Randomized Linear Algebra, over at <https://arxiv.org/pdf/1608.04481>.

The Wikipedia page on McDoarmid's inequality, over at <https://en.wikipedia.org/wiki/McDiarmid%27s_inequality>.