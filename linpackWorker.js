onmessage = function(event) {
    "use strict";
    
    var ret = event.data;
    var asize = ret.arraySize;
    var ln = new Linpack();

    ln.setArsize(asize);
    var text = ln.runLinpack();
    
    postMessage(text);
}

function Linpack() {
    "use strict";
    var DoubleLength = 64;
    var IntegerLength = 32;
    var arsize = 200;
    var result;
}

Linpack.prototype.setArsize = function(a) {
    "use strict";
    a /= 2;
    a *= 2;
    this.arsize = a;
    if (this.arsize < 10) {alert("Too Small");}
};    

Linpack.prototype.getArsize = function() {
    return this.arsize;
};

Linpack.prototype.runLinpack = function()
{
    "use strict";
    var arsize2d, memreq,nreps;   
    this.arsize/=2;
    this.arsize*=2;

    arsize2d = this.arsize * this.arsize;
    memreq = arsize2d * this.DoubleLength + this.arsize * this.DoubleLength + this.arsize * this.IntegerLength;
//    alert("Memory required: " + (memreq+512) * 2* 2* 10);
    nreps=1;
    while (this.linpackSub(nreps, this.arsize)/1000 <5 )
    {
        nreps *= 2;
    }
    return(this.result);
};

Linpack.prototype.linpackSub = function(nreps, arsize)
{
    "use strict";
    var a = [];
    var b = [];
    var norma, kflops, ops;
    var ipvt = [];
    var n,info,lda;
    var i, arsize2d;
    var start, t1;
    var tdgefa, tdgesl, totalt, toverhead;
    norma = 0;
    info = 0;
    this.result = "";

    lda = arsize;
    n = arsize/2;
    arsize2d =  arsize * arsize;
    ops=((2.0*n*n*n)/3.0+2.0*n*n);
    a= [arsize2d];
    b= [arsize]; 
    ipvt= [arsize];
    tdgesl= 0;
    tdgefa= 0;
    start= new Date().getTime();
    for (i=0;i<nreps;i++)
        {
        norma = this.matgen(a,lda,n,b,norma);
        t1 = new Date().getTime();
        info = this.dgefa(a,lda,n,ipvt, info,1);
        tdgefa+=(new Date().getTime()-t1);
        t1 = new Date().getTime();
        this.dgesl(a,lda,n,ipvt,b,0,1);
        tdgesl += (new Date().getTime()-t1);
        }
    for (i=0;i<nreps;i++)
        {
        norma = this.matgen(a,lda,n,b,norma);
        t1 = new Date().getTime();
        info = this.dgefa(a,lda,n,ipvt, info,0);
        tdgefa += new Date().getTime()-t1;
        t1 = new Date().getTime();
        this.dgesl(a,lda,n,ipvt,b,0,0);
        tdgesl += new Date().getTime()-t1;
        }
    totalt= new Date().getTime()-start;
    if (totalt/1000<0.5 || (tdgefa+tdgesl)/1000<0.2) {
        return(0);
    }
    kflops=2*nreps*ops/(1000*(tdgefa+tdgesl)/1000);
    toverhead=totalt-tdgefa-tdgesl;
    if (tdgefa/1000<0)
        {tdgefa=0;}
    if (tdgesl/1000<0)
        {tdgesl=0;}
    if (toverhead/1000<0)
        {toverhead=0;}
    // for (int j = 0; j < b.Length; j++)
    // {
    //    Console.WriteLine(b[j]);
    // }
    this.result +=
        nreps + " " + ((totalt/1000)) + " " + (100 * tdgefa / totalt) +
        " " + (100 * tdgesl / totalt) + " " + (100 * toverhead / totalt) +
        " " +(kflops);
    return(totalt);
};

/*
** For matgen,
** We would like to declare a[][lda], but c does not allow it.  In this
** function, references to a[i][j] are written a[lda*i+j].
*/
Linpack.prototype.matgen = function(a, lda, n, b, norma)
{
    "use strict";
    var init, i, j;
    init = 1325;
    norma = 0.0;
    for (j = 0; j < n; j++)
    {
        for (i = 0; i < n; i++)
        {
            init = 3125 * init % 65536;
            a[lda*j+ i] = (init - 32768.0) / 16384.0;
            norma = (a[lda*j+ i] > norma) ? a[lda*j+ i] : norma;
        }
    }
    for (i = 0; i < n; i++)
    {
        b[i] = 0.0;
    }
    for (j = 0; j < n; j++)
    {
        for (i = 0; i < n; i++) {b[i] = b[i] + a[lda*j+ i];}
    }
    return norma;
};

/*
**
** DGEFA benchmark
**
** We would like to declare a[][lda], but c does not allow it.  In this
** function, references to a[i][j] are written a[lda*i+j].
**
**   dgefa factors a double precision matrix by gaussian elimination.
**
**   dgefa is usually called by dgeco, but it can be called
**   directly with a saving in time if  rcond  is not needed.
**   (time for dgeco) = (1 + 9/n)*(time for dgefa) .
**
**   on entry
**
**      a       REAL precision[n][lda]
**              the matrix to be factored.
**
**      lda     integer
**              the leading dimension of the array  a .
**
**      n       integer
**              the order of the matrix  a .
**
**   on return
**
**      a       an upper triangular matrix and the multipliers
**              which were used to obtain it.
**              the factorization can be written  a = l*u  where
**              l  is a product of permutation and unit lower
**              triangular matrices and  u  is upper triangular.
**
**      ipvt    integer[n]
**              an integer vector of pivot indices.
**
**      info    integer
**              = 0  normal value.
**              = k  if  u[k][k] .eq. 0.0 .  this is not an error
**                   condition for this subroutine, but it does
**                   indicate that dgesl or dgedi will divide by zero
**                   if called.  use  rcond  in dgeco for a reliable
**                   indication of singularity.
**
**   linpack. this version dated 08/14/78 .
**   cleve moler, university of New Mexico, argonne national lab.
**
**   functions
**
**   blas daxpy,dscal,idamax
**
*/

Linpack.prototype.dgefa = function(a, lda, n, ipvt, info, roll)
{  
    "use strict";
    var t;
    var j, k, kp1, l,nm1;

    /* gaussian elimination with partial pivoting */
    if (roll !== 0)
    {  
        info = 0;
        nm1 = n - 1;
        if (nm1 >=  0)
        {  
            for (k = 0; k < nm1; k++)
            {  
                kp1 = k + 1;

                /* find l = pivot index */
                l = this.idamax(n-k,a,lda*k+k,1) + k;
                ipvt[k] = l;

                /* zero pivot implies this column already
                    triangularized */
                if (a[lda*k+l] !== 0)
                {  
                    /* interchange if necessary */
                    if (l !== k)
                    {  
                        t = a[lda*k+l];
                        a[lda*k+l] = a[lda*k+k];
                        a[lda*k+k] = t;
                    }

                    /* compute multipliers */

                    t = -1/a[lda*k+k];
                    this.dscal_r(n-(k+1),t,a,lda*k+k+1,1);

                    /* row elimination with column indexing */

                    for (j = kp1; j < n; j++)
                    {  
                        t = a[lda*j+l];
                        if (l !== k)
                        {  
                            a[lda*j+l] = a[lda*j+k];
                            a[lda*j+k] = t;
                        }
                        this.daxpy_r(n-(k+1),t,a, lda*k+k+1,1,a,lda*j+k+1,1);
                    }
                }
                else
                    {info = k;}
            }
        }
        ipvt[n-1] = n-1;
        if (a[lda*(n-1)+(n-1)] === 0) {info = n-1;}
    }
    else
    {
        info = 0;
        nm1 = n - 1;
        if (nm1 >=  0)
        {
            for (k = 0; k < nm1; k++)
            {
                kp1 = k + 1;

                /* find l = pivot index */
                l = this.idamax(n-k,a,lda*k+k,1) + k;
                ipvt[k] = l;

                /* zero pivot implies this column already
                    triangularized */
                if (a[lda*k+l] !== 0)
                {
                    /* interchange if necessary */
                    if (l !== k)
                    {
                        t = a[lda*k+l];
                        a[lda*k+l] = a[lda*k+k];
                        a[lda*k+k] = t;
                    }

                    /* compute multipliers */
                    t = -1/a[lda*k+k];
                    this.dscal_ur(n-(k+1),t,a,lda*k+k+1,1);

                    /* row elimination with column indexing */
                    for (j = kp1; j < n; j++)
                    {
                        t = a[lda*j+l];
                        if (l !== k)
                        {  
                            a[lda*j+l] = a[lda*j+k];
                            a[lda*j+k] = t;
                        }
                        this.daxpy_ur(n-(k+1),t,a,lda*k+k+1,1,a,lda*j+k+1,1);
                    }
                }
                else
                    {info = k;}
            }
        }
        ipvt[n-1] = n-1;
              if (a[lda*(n-1)+(n-1)] === 0) {info = n-1;}
    }
    return info;
};

/*
**
** DGESL benchmark
**
** We would like to declare a[][lda], but c does not allow it.  In this
** function, references to a[i][j] are written a[lda*i+j].
**
**   dgesl solves the double precision system
**   a * x = b  or  trans(a) * x = b
**   using the factors computed by dgeco or dgefa.
**
**   on entry
**
**      a       double precision[n][lda]
**              the output from dgeco or dgefa.
**
**      lda     integer
**              the leading dimension of the array  a .
**
**      n       integer
**              the order of the matrix  a .
**
**      ipvt    integer[n]
**              the pivot vector from dgeco or dgefa.
**
**      b       double precision[n]
**              the right hand side vector.
**
**      job     integer
**              = 0         to solve  a*x = b ,
**              = nonzero   to solve  trans(a)*x = b  where
**                          trans(a)  is the transpose.
**
**  on return
**
**      b       the solution vector  x .
**
**   error condition
**
**      a division by zero will occur if the input factor contains a
**      zero on the diagonal.  technically this indicates singularity
**      but it is often caused by improper arguments or improper
**      setting of lda .  it will not occur if the subroutines are
**      called correctly and if dgeco has set rcond .gt. 0.0
**      or dgefa has set info .eq. 0 .
**
**   to compute  inverse(a) * c  where  c  is a matrix
**   with  p  columns
**         dgeco(a,lda,n,ipvt,rcond,z)
**         if (!rcond is too small){
**              for (j=0,j<p,j++)
**                      dgesl(a,lda,n,ipvt,c[j][0],0);
**         }
**
**   linpack. this version dated 08/14/78 .
**   cleve moler, university of new mexico, argonne national lab.
**
**   functions
**
**   blas daxpy,ddot
*/
Linpack.prototype.dgesl = function(a, lda, n, ipvt, b, job, roll)
{  
    "use strict";
    var t;
    var k, kb, l, nm1;

    if (roll !== 0)
    {  
        nm1 = n - 1;
        if (job === 0)
        {  
            /* job = 0 , solve  a * x = b   */
            /* first solve  l*y = b         */
            if (nm1 >= 1)
            {  
                for (k = 0; k < nm1; k++)
                {  
                    l = ipvt[k];
                    t = b[l];
                    if (l !== k)
                    {  
                        b[l] = b[k];
                        b[k] = t;
                    }
                    this.daxpy_r(n - (k + 1), t, a,lda * k + k + 1, 1, b,k + 1, 1);
                }
            }
            /* now solve  u*x = y */
            for (kb = 0; kb < n; kb++)
            {  
                k = n - (kb + 1);
                b[k] = b[k] / a[lda * k + k];
                t = -b[k];
                this.daxpy_r(k, t, a, lda*k+0, 1, b, 0, 1);
            }
        }
        else
        {  
            /* job = nonzero, solve  trans(a) * x = b  */
            /* first solve  trans(u)*y = b             */
            for (k = 0; k < n; k++)
            {  
                t = this.ddot_r(k, a, lda * k + 0, 1, b,0, 1);
                b[k] = (b[k] - t) / a[lda * k + k];
            }
            /* now solve trans(l)*x = y     */

            if (nm1 >= 1)
            {  
                for (kb = 1; kb < nm1; kb++)
                {  
                    k = n - (kb + 1);
                    b[k] = b[k] + this.ddot_r(n - (k + 1), a, lda * k + k + 1, 1, b, k + 1, 1);
                    l = ipvt[k];
                    if (l !== k)
                    {  
                        t = b[l];
                        b[l] = b[k];
                        b[k] = t;
                    }
                }
            }
        }
    }
    else
    {  
        nm1 = n - 1;
        if (job === 0)
        {  
            /* job = 0 , solve  a * x = b   */
            /* first solve  l*y = b         */
            if (nm1 >= 1)
            {  
                for (k = 0; k < nm1; k++)
                {  
                    l = ipvt[k];
                    t = b[l];
                    if (l !== k)
                    {  
                        b[l] = b[k];
                        b[k] = t;
                    }
                    this.daxpy_ur(n - (k + 1), t, a, lda * k + k + 1, 1, b, k + 1, 1);
                }
            }
            /* now solve  u*x = y */
            for (kb = 0; kb < n; kb++)
            {  
                k = n - (kb + 1);
                b[k] = b[k] / a[lda * k + k];
                t = -b[k];
                this.daxpy_ur(k, t, a, lda * k + 0, 1, b, 0, 1);
            }
        }
        else
        {  
            /* job = nonzero, solve  trans(a) * x = b  */
            /* first solve  trans(u)*y = b             */
            for (k = 0; k < n; k++)
            {  
                t = this.ddot_ur(k, a, lda * k + 0, 1, b, 0, 1);
                b[k] = (b[k] - t) / a[lda * k + k];
            }

            /* now solve trans(l)*x = y     */
            if (nm1 >= 1)
            {  
                for (kb = 1; kb < nm1; kb++)
                {  
                    k = n - (kb + 1);
                    b[k] = b[k] + this.ddot_ur(n - (k + 1), a, lda * k + k + 1, 1, b,k + 1, 1);
                    l = ipvt[k];
                    if (l !== k)
                    {  
                        t = b[l];
                        b[l] = b[k];
                        b[k] = t;
                    }
                }
            }
        }
    }
};

/*
** Constant times a vector plus a vector.
** Jack Dongarra, linpack, 3/11/78.
** ROLLED version
*/
Linpack.prototype.daxpy_r = function(n, da, dx, offset_dx, incx, dy, offset_dy, incy)
{  
    "use strict";
    var i, ix, iy;
    var of_dx = offset_dx;
    var of_dy = offset_dy;

    if (n <= 0) {return;}
    if (da === 0) {return;}

    if (incx !== 1 || incy !== 1)
    {  
        /* code for unequal increments or equal increments != 1 */
        ix = 1;
        iy = 1;
        if (incx < 0) {ix = (-n + 1) * incx + 1;}
        if (incy < 0) {iy = (-n + 1) * incy + 1;}
        for (i = 0; i < n; i++)
        {  
            dy[iy+of_dy] = dy[iy+of_dy] + da * dx[ix+of_dx];
            ix = ix + incx;
            iy = iy + incy;
        }
        return;
    }

    /* code for both increments equal to 1 */
    for (i = 0; i < n; i++) {dy[i+of_dy] = dy[i+of_dy] + da * dx[i+of_dx];}
};

/*
** Forms the dot product of two vectors.
** Jack Dongarra, linpack, 3/11/78.
** ROLLED version
*/
Linpack.prototype.ddot_r = function(n, dx, offset_dx, incx, dy, offset_dy, incy)
{  
    "use strict";
    var dtemp;
    var i, ix, iy;
    var of_dx = offset_dx;
    var of_dy = offset_dy;

    dtemp = 0;

    if (n <= 0) {return (0);}

    if (incx !== 1 || incy !== 1)
    {  
        /* code for unequal increments or equal increments != 1 */
        ix = 0;
        iy = 0;
        if (incx < 0) {ix = (-n + 1) * incx;}
        if (incy < 0) {iy = (-n + 1) * incy;}
        for (i = 0; i < n; i++)
        {  
            dtemp = dtemp + dx[ix+of_dx] * dy[iy+of_dy];
            ix = ix + incx;
            iy = iy + incy;
        }
        return (dtemp);
    }

    /* code for both increments equal to 1 */
    for (i = 0; i < n; i++) {dtemp = dtemp + dx[i+of_dx] * dy[i+of_dy];}
    return (dtemp);
};

/*
** Scales a vector by a constant.
** Jack Dongarra, linpack, 3/11/78.
** ROLLED version
*/
Linpack.prototype.dscal_r = function(n, da, dx, offset_dx, incx)
{
    "use strict";
    var i, nincx;
    var of_dx = offset_dx;

    if (n <= 0) {return;}
    if (incx !== 1)
    {
        /* code for increment not equal to 1 */
        nincx = n * incx;
        for (i = 0; i < nincx; i = i + incx) {dx[i+of_dx] = da * dx[i+of_dx];}
        return;
    }

    /* code for increment equal to 1 */
    for (i = 0; i < n; i++) {dx[i+of_dx] = da * dx[i+of_dx];}
};

/*
** constant times a vector plus a vector.
** Jack Dongarra, linpack, 3/11/78.
** UNROLLED version
*/
Linpack.prototype.daxpy_ur = function(n, da, dx, offset_dx, incx, dy, offset_dy, incy)
{  
    "use strict";
    var i, ix, iy, m;
    var of_dx = offset_dx;
    var of_dy = offset_dy;
    if (n <= 0) {return;}
    if (da === 0) {return;}

    if (incx !== 1 || incy !== 1)
    {  
        /* code for unequal increments or equal increments != 1 */
        ix = 1;
        iy = 1;
        if (incx < 0) {ix = (-n + 1) * incx + 1;}
        if (incy < 0) {iy = (-n + 1) * incy + 1;}
        for (i = 0; i < n; i++)
        {  
            dy[iy+of_dx] = dy[iy+of_dy] + da * dx[ix+of_dx];
            ix = ix + incx;
            iy = iy + incy;
        }
        return;
    }

    /* code for both increments equal to 1 */
    m = n % 4;
    if (m !== 0)
    {  
        for (i = 0; i < m; i++) {dy[i+of_dy] = dy[i+of_dy] + da * dx[i+of_dx];}
        if (n < 4) {return;}
    }
    for (i = m; i < n; i = i + 4)
    {  
        var i_dx = i + of_dx;
        var i_dy = i + of_dy;
        dy[i_dy] = dy[i_dy] + da * dx[i_dx];
        dy[i_dy + 1] = dy[i_dy + 1] + da * dx[i_dx + 1];
        dy[i_dy + 2] = dy[i_dy + 2] + da * dx[i_dx + 2];
        dy[i_dy + 3] = dy[i_dy + 3] + da * dx[i_dx + 3];
    }
};

/*
** Forms the dot product of two vectors.
** Jack Dongarra, linpack, 3/11/78.
** UNROLLED version
*/
Linpack.prototype.ddot_ur = function(n, dx, offset_dx, incx, dy, offset_dy, incy)
{  
    "use strict";
    var dtemp;
    var of_dx = offset_dx;
    var of_dy = offset_dy;
    var i, ix, iy, m;

    dtemp = 0;

    if (n <= 0) {return (0);}

    if (incx !== 1 || incy !== 1)
    {  
        /* code for unequal increments or equal increments != 1 */
        ix = 0;
        iy = 0;
        if (incx < 0) {ix = (-n + 1) * incx;}
        if (incy < 0) {iy = (-n + 1) * incy;}
        for (i = 0; i < n; i++)
        {  
            dtemp = dtemp + dx[ix+of_dx] * dy[iy+of_dy];
            ix = ix + incx;
            iy = iy + incy;
        }
        return (dtemp);
    }

    /* code for both increments equal to 1 */
    m = n % 5;
    if (m !== 0)
    {  
        for (i = 0; i < m; i++) {dtemp = dtemp + dx[i+of_dx] * dy[i+of_dy];}
        if (n < 5) {return (dtemp);}
    }
    for (i = m; i < n; i = i + 5)
    {  
        var i_dx = i + of_dx;
        var i_dy = i + of_dy;
        dtemp = dtemp + dx[i_dx] * dy[i_dy] +
        dx[i_dx + 1] * dy[i_dy + 1] + dx[i_dx + 2] * dy[i_dy + 2] +
        dx[i_dx + 3] * dy[i_dy + 3] + dx[i_dx + 4] * dy[i_dy + 4];
    }
    return (dtemp);
};

/*
** Scales a vector by a constant.
** Jack Dongarra, linpack, 3/11/78.
** UNROLLED version
*/
Linpack.prototype.dscal_ur = function(n, da, dx, offset_dx, incx)
{  
    "use strict";
    var i, m, nincx;
    var of_dx = offset_dx;
    if (n <= 0) {return;}
    if (incx !== 1)
    {  
        /* code for increment not equal to 1 */
        nincx = n * incx;
        for (i = 0; i < nincx; i = i + incx) {dx[i+of_dx] = da * dx[i+of_dx];}
        return;
    }

    /* code for increment equal to 1 */
    m = n % 5;
    if (m !== 0)
    {  
        for (i = 0; i < m; i++) {dx[i+of_dx] = da * dx[i+of_dx];}
        if (n < 5) {return;}
    }
    for (i = m; i < n; i = i + 5)
    {  
        var i_tmp = i + of_dx;
        dx[i_tmp] = da * dx[i_tmp];
        dx[i_tmp + 1] = da * dx[i_tmp + 1];
        dx[i_tmp + 2] = da * dx[i_tmp + 2];
        dx[i_tmp + 3] = da * dx[i_tmp + 3];
        dx[i_tmp + 4] = da * dx[i_tmp + 4];
    }
};

/*
** Finds the index of element having max. absolute value.
** Jack Dongarra, linpack, 3/11/78.
*/
Linpack.prototype.idamax = function(n, dx, offset_dx, incx)
{  
    "use strict";
    var dmax;
    var of_dx = offset_dx;
    var ix, itemp, i;
    itemp = 0;
    if (n < 1)
            {return (-1);}
    if (n === 1)
            {return (0);}
    if (incx !== 1)
    {
        /* code for increment not equal to 1 */
        ix = 1;
        dmax = Math.abs((dx[0+of_dx]));
        ix = ix + incx;
        for (i = 1; i < n; i++)
        {
            if (Math.abs(dx[ix+of_dx]) > dmax)
            {
                itemp = i;
                dmax = Math.abs(dx[ix+of_dx]);
            }
            ix = ix + incx;
        }
    }
    else
    {
        /* code for increment equal to 1 */
        itemp = 0;
        dmax = Math.abs(dx[0+of_dx]);
        for (i = 1; i < n; i++)
        {
            if (Math.abs(dx[i+of_dx]) > dmax)
            {
                itemp = i;
                dmax = Math.abs(dx[i+of_dx]);
            }
        }
    }
    return (itemp);
};

