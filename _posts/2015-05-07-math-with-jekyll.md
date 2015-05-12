---
layout: post
title:  "Math with Jekyll"
date:   2015-05-07
categories: jekyll markdown kramdown math octave
---


The goal
========
Hopefully, by using [Jekyll](http://jekyllrb.com), [kramdown](http://kramdown.gettalong.org), and [MathJax](http://www.mathjax.org), I should be able to get pretty-looking math to show up here.


The attempt
===========
In the standard Jekyll directory structure, edit the file **_include/head.html** file. Add the following inside the `<head> ... </head>` tags:

`<script type="text/javascript" src="//cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"></script>`

That's it. Now, whenever you want to write some $$\LaTeX$$ style formula, just enclose it in `$$ ... $$`.

The test
========
We initialise a scalar, e.g. $$ \alpha = \frac{1}{2} $$, and then use it as a relaxation parameter for [Newton's method](https://en.wikipedia.org/wiki/Newton's_method) to search for a root:

$$
\begin{align*}
	x_{i+1} = x_{i} - \alpha \frac{f(x_i)}{f^{\prime}(x_i)}.
\end{align*}
$$

If we start with $$ x_0 = 1 $$ and $$ f(x) = 4 - x^2 $$, the derivative is $$ f^{\prime}(x) = -2x $$, we can get successive values for $$ x $$ that approach a nearby root of $$f$$.

The first few points (rounded to three decimal places):

|  i |   x   |   y   |
|----+-------+-------|
|  0 | 1.000 | 3.000 |
|  1 | 1.750 | 0.938 |
|  2 | 1.884 | 0.451 |
|  3 | 1.944 | 0.222 |
|  4 | 1.972 | 0.110 |
|  5 | 1.986 | 0.055 |
|  6 | 1.993 | 0.027 |
|  7 | 1.997 | 0.014 |
|  8 | 1.998 | 0.007 |
|  9 | 1.999 | 0.003 |
| 10 | 2.000 | 0.002 |
| 11 | 2.000 | 0.001 |



If instead we wanted to find a local optima, we might instead try to find the root of the first derivative. In this case, we will need the second derivative of our function. This is obviously $$ f^{\prime\prime}(x) = -2 $$. Our iterative step is now

$$
x_{i+1} = x_{i} - \alpha \frac{f^{\prime}(x)}{f^{\prime\prime}(x)} = x_{i} - \frac{x_{i}}{2} = \frac{x_{i}}{2}.
$$

As we iterate, we now get:

|  i |   x   |   y   |
|----+-------+-------|
|  0 | 1.000 | 3.000 |
|  1 | 0.500 | 3.750 |
|  2 | 0.250 | 3.938 |
|  3 | 0.125 | 3.984 |
|  4 | 0.063 | 3.996 |
|  5 | 0.031 | 3.999 |


The code
========
The above values were generated with the following Octave code:

~~~matlab
function y = f(x), y = 4 - x**2; end;
function dy = fp(x), dy = -2*x; end;
function d2y = fpp(x), d2y = -2; end;
# Find a root.
x=1; i=0;
while (abs(f(x)) > 0.001)
  i = i + 1;
  x = x - 0.5 * f(x) / fp(x);
  printf('%i %f %f\n', i, x, f(x));
end
# Find a minima/maxima.
x=1; i=0;
while (abs(4-f(x)) > 0.001)
  i = i + 1;
  x = x - 0.5 * fp(x) / fpp(x);
  printf('%i %f %f\n', i, x, f(x));
end
~~~

Note that the omniscient termination criterion locating the maxima was used to avoid the remaining points all sitting at $$ y = 4.000 $$. If the while loop instead read `while (abs(fp(x)) > 0.001)`, then it would have taken 11 iterations to converge, halving the value of $$ x $$ at each iteration as expected.

Conclusion
==========
Human-readable text files, pretty looking math, easy.


