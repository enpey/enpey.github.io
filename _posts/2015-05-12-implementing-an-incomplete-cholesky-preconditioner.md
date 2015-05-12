---
layout: post
title: "Implementing an incomplete Cholesky preconditioner"
date: 2015-05-12
categories: preconditioning Cholesky factorisation
---


Incomplete Cholesky factorisations form efficient preconditioners for solving well-conditioned symmetric positive definite (SPD) linear systems with solvers like the [Conjugate Gradient method](https://en.wikipedia.org/wiki/Conjugate_gradient_method). These preconditioners are fairly robust and provide an alternative to full factorisation which may lead to a reduction in the time required to solve the given linear system of equations. In the following, an efficient incomplete Cholesky factorisation for sparse linear systems is outlined. The method is a left-looking factorisation based on the description given by [Timothy Davis](http://faculty.cse.tamu.edu/davis/welcome.html) in *Direct Methods for Sparse Linear Systems*.


## Overview ##
The factorisation is *incomplete* as it drops entries during the factorisation. Common dropping rules rely on limiting the number of entries per row/column or dropping entries that are deemed to be small in magnitude based on some threshold value. (There are level-of-fill based methods, too, but these are not considered here.) In the method described here, both fill and threshold dropping are used. The symmetry of the coefficient matrix is exploited by only requiring the diagonal and strict lower triangle.


### Input ###
The following input is required:
  
  * scaling vector, `scal`, that has n entries that will symetrically scale $$\mathbf{A}$$, i.e. if $$\mathbf{S} = \text{diag}(\text{scal})$$, then we factorise $$\mathbf{A}_{s} = \mathbf{SAS}$$;
  * diagonal shift parameters, `alpha` and `beta`, that modify the matrix diagonal: $$\mathbf{A}_{ss} = \mathbf{A}_{s} + \alpha \text{diag}(\mathbf{A}_{s}) + \beta \mathbf{I}$$;
  * the coefficient matrix, $$\mathbf{A}$$, in compressed sparse column (CSC) format (following Davis, `Ap`, `Ai`, and `Ax`);
  * the size of $$A$$, `m` $$\times$$ `n`;
  * the drop threshold, `tau`; and
  * the maximum number of entries per column, `Lfill`.


### Output ###
The following output is provided:

  * the incomplete factor, $$\mathbf{L}$$, in CSC format (`Lp`, `Li`, `Lx`).


### Local variables ###
We will need a range of local variables to keep track of various things. The major variables are:

  * `head`: an integer n-vector that will point to the beginning of one of n linked lists (for the update);
  * `next`: an integer n-vector that points to the next component in one of the n n linked lists;
  * `Lpos`: an integer n-vector that points to the next location in `Lx` to assist the update procedure;
  * `flag`: an integer n-vector that will flag which rows already have an entry in the column being computed (this uses `next`, but is different here for clarity);
  * `Vx`: a real n-vector that will hold the non-zero values of the column being computed;
  * `Vi`: an integer n-vector that will hold the row indices of the non-zero entries in `Vx`;
  * `d`: a real n-vector that will hold the values of the diagonal entries (this will enable easy comparison against any diagonal value in our dropping scheme, as well as allowing us to halt as soon as any diagonal entry becomes less than or equal to zero);


## Method description ##
The method is broken into five sections:

  1. Initialisation
  2. Gather and scaling
  3. Update
  4. Dropping
  5. Divide and store


The last three sections are contained within the main loop and will be executed over and over again.


### Initialisation ###
The only vectors that need initialising are `d`, `head`, and `flag`. We set `head` to $$-1$$, `flag` to $$0$$ and `d` to the diagaonl entries of $$\mathbf{A}$$, scaling them as necessary.

~~~fortran
 do j = 1, n
    flag(j) = 0
    head(j) = -1
    do ii = Ap(j), Ap(j+1)-1
       if (Ai(ii) == j) then
          d(j) = (alpha + one) * a(ii) * scal(j)**2 + beta
	  exit
       end if
    end do
 end do

 Lp(1) = 1
 Lnz = 0    ! Number of entries in L.
 Lextra = 0 ! Allowed fill entries that are unused.
~~~

Following initialisation, the main loop is started and loops from `k = 1` to `n`.


### Gather and scaling ###
The first thing we do for each `k` is gather the $$k$$'th column of $$\mathbf{A}$$ into our vectors `Vx`, `Vi`, and `flag`. We also take the opportunity to scale the column here. We use a counter to keep track of how many non-zero entries we have in this column, `colnz`.

~~~fortran
 do ii = Ap(k), Ap(k+1)-1
    i = Ai(ii)
    if (i > k) then
       colnz = colnz + 1
       Vx(i) = scal(i) * Ax(ii) * scal(k)
       Vi(colnz) = i
       flag(i) = k
    end if
 end do
~~~


### Update ###
The update forms the bulk of the method. All we are doing is walking along the $$k$$'th row of $$\mathbf{L}$$ and, wherever there is a non-zero entry, $$l_{kj}$$, subtracting $$l_{kj} \times \mathbf{L}_{j}$$ from the column we are currently computing (where $$\mathbf{L}_{j}$$ is the $$j$$'th column of $$\mathbf{L}$$). This is where the linked list is used.

We use the linked list to walk along the $$k$$'th row without doing any messy searching through the CSC data structures. We start with the non-zero entry in the $$k$$'th row in column `j = head(k)`. The next entry is in the column `next(j)`, and so on, until `next(j) == -1`, indicating the end of the row. As we move along the row, we also need to update the linked lists required for later updates after finishing with each $$l_{kj}$$. We do this incrementing `Lpos(j)` and checking whether there any entries remaining in the $$j$$'th column. If there is at least one non-zero entry remaining in the column, let us call its row index `ancestor`. We then push the $$j$$'th column onto the head of the `ancestor`'th linked list.

The update procedure happens before updating the linked list. The update simply walks through the $$j$$'th column of $$\mathbf{L}$$ below row $$k$$, and subtracts $$l_{kj}$$ times the entries from our current column. For each subtraction though, we need to check if we already have an entry or not in our column. This is easily checked by considering whether `flag(i) == k` or not ($$i$$ is the row index of the entry we have found in column $$j$$). If it is, we already have an entry in position $$(i,k)$$ and just need to update it. If we do not have an entry, then we must increment `colnz`, set `flag(i) = k`, and add the entry in our `Vx` an `Vi` vectors.

All together, we get the following.

~~~fortran
 j = head(k)
 do while (j /= -1)
    lkj = Lx(Lpos(j)) ! Get L(k,j).
    do ii = Lpos(j)+1, Lp(j+1)-1 ! Walk through column j below row k.
       i = Li(ii)
       if (flag(i) == k) then ! Already have an entry in L(i,k).
	  Vx(i) = Vx(i) - lkj * Lx(ii)
       else ! We have fill-in.
	  colnz = colnz + 1
	  flag(i) = k
	  Vx(i) = -lkj * Lx(ii)
	  Vi(colnz) = i
       end if
    end do

    ! Now, update the linked list.
    jnext = next(j)
    Lpos(j) = Lpos(j) + 1
    if (Lpos(j) < Lp(j+1)) then
       ancestor = Li(Lpos(j))
       next(j) = head(ancestor)
       head(ancestor) = j
    end if
    j = jnext
 end do
~~~


### Dropping ###
Now we decide which entries in our column (held in `Vx` and `Vi`) we should keep. Note that we cannot necessarily assume there is enough space allocated in our $$\mathbf{L}$$ structures to hold the entries we keep, unless we are using `Lfill` to restrict the number of entries in each column (in which case the $$\mathbf{L}$$ structures should have been pre-allocated to the correct size). We allow our method to add as many entries as it needs to if `Lfill < 0` by ensuring that if we were able to store all `colnz` entries in bour $$\mathbf{L}$$ structures (if there is not enough space in `Li` and `Lx`, we reallocate them with more space).

The first round of dropping uses the threshold `tau`. We compare each entry in our column against `tau` $$\times$$ `dsqrt(d(k))`. If it is smaller than or equal to this threshold, we drop it. If it is greater than it, we store in $$\mathbf{L}$$.

~~~fortran
 if (tau > zero) then
    drop = tau * dsqrt(d(k))
    colnz = 0
    do ii = 1, colnz0
       if (abs(Vx(ii)) > drop) then
	  colnz = colnz + 1
	  Vi(colnz) = Vi(ii)
       end if
    end do
 end if
~~~

If we are using fill control, we need to make sure that there are no more than `Lfill` entries in our column. We must consider the following cases:

  1. `Lfill < 0`: no fill control so we are keeping all entries that passed the threshold test;
  2. `colnz <= Lfill`: we can keep all the entries that passed the threshold test; or
  3. `colnz > Lfill`: we must find the `Lfill` largest entries in the column and discard the rest.

It is only the third case where we need to do anything. We use a *quicksplit*-style method (described by [Yousef Saad](http://www-users.cs.umn.edu/~saad/) in *Iterative Methods for Sparse Linear Systems*) to extract the largest `Lfill` entries in the column (this method is based on the *quicksort* algorithm and avoids sorting the entries and then taking the desired entries from the end). Note that in the second case we end up with unused space. Obviously, it will be beneficial if we can use this space later on. We use `Lextra` to keep track of any extra space that arises from this case and allow all extra space to be used beginning with the next column.

The final step is to sort the column's entries into ascending row order. This is necessary to ensure that our linked list update procedure works. Note that the linked list is not the only, nor even necessarily the fastest, way to accomplish this (efficient alternatives include graph-based methods and binary search trees). For this, we use quicksort.

~~~fortran
 if (Lfill < 0) then
    allowed = colnz
 else
    allowed = Lfill + Lextra ! In case there is some unused space from an earlier column.
 end if
 if (colnz > allowed) then
    call quicksplit_lrg_ind(colnz, Vi, Vx, allowed)
    colnz = allowed
 end if
 Lextra = allowed - colnz ! Keep track of any unused space.
 call quicksort(colnz, Vi) ! Sort indices into ascending row order.
~~~


### Divide and store ###
All that is left is to divide the column by the square root of the diagonal and store the column in our $$\mathbf{L}$$ structures. We also update the diagonal entries by subtracting the square of each entry in our column and checking that those diagonal entries modified stay positive, and then we push the first off-diagonal entry onto the associated linked list.

~~~fortran
 lkk = d(k)
 if (lkk <= zero) then
    ! Indicate failure and return.
    info(1) = k
    return
 end if

 lkk = dsqrt(lkk)
 Lnz = Lnz + 1
 Lx(Lnz) = lkk
 Li(Lnz) = k
 if (colnz /= 0) then
    lkk = one / lkk
    do ii = 1, colnz
       Lnz = Lnz + 1
       i = Vi(ii)
       Li(Lnz) = i
       Lx(Lnz) = Vx(i) * lkk
       d(i) = d(i) - Lx(Lnz)**2
       if (d(i) <= zero) then
	  info(1) = i
	  return
       end if
    end do

    ! Add to linked list for parent.
    Lpos(k) = Lp(k) + 1
    ancestor = Li(Lpos(k))
    next(k) = head(ancestor)
    head(ancestor) = k
 end if
 head(k) = -1
~~~


## Summary ##
The above described an efficient implementation of an incomplete Cholesky factorisation. It allows the user to input diagonal shifts to increase the diagonal dominance of the matrix and symmetric scaling. The dropping criteria consists of a threshold value and a fill control, both of which may be set such that they have no influence on the factorisation.