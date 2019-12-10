mpz_t euclidian(mpz_t a) /* Extended Euclidian algorithm */
{
	mpz_t u, v, g, u1, v1, g1, q, t1, t1_, t2, t2_, t3, t3_, rop;

	mpz_init_set_ui(u,1);
	mpz_init_set_ui(v,0);
	mpz_init_set(g, a);
	mpz_init_set_ui(u1,0);
	mpz_init_set_ui(v1,1);
	mpz_init_set(g1, ec.p);
	mpz_inits(q, t1, t1_, t2, t2_, rop);

	//~ while (g1 != 0)
	while(mpz_cmp_ui(g1, 0) != 0)
	{
		//~ q = g / g1;
		mpz_cdiv_q(q, g, g1); 
		
		//~ t1 = u - q * u1;
		mpz_mul(t1_, q, u1);
		mpz_sub(t1, u, t1_);
		
		//~ t2 = v - q * v1;
		mpz_mul(t2_, q, v1);
		mpz_sub(t2, v, t2_);

		//~ t3 = g - q * g1;
		mpz_mul(t3_, q, g1);
		mpz_sub(t3, g, t3_);

		mpz_set(u, u1);
		mpz_set(v, v1);
		mpz_set(g, g1);
		mpz_set(u1, t1);
		mpz_set(v1, t2);
		mpz_set(g1, t3);
	}
	
	mpz_mod(rop, u, ec.n);
	mpz_clears(u, v, g, u1, v1, g1, q, t1, t1_, t2, t2_);
	return rop;
}

mpz_t mod(mpz_t a)
{
	mpz_t rop, tmp;
	mpz_init(rop);	
	mpz_init_set(tmp, a);
	if(mpz_cmp(a, 0) < 0)
		mpz_add(tmp, a, ec.n);
	mpz_mod(rop, tmp, ec.n);
	mpz_clear(tmp);
	return rop;
}

void addInv(struct point *q) /* Additive inverse*/
{
	//~ P - Q = P + (-Q)
	mpz_t rop;
	mpz_init(rop);
	mpz_mul_si(rop, q->y, -1);
	
	q->y = mod(q->y * -1); 
}

mpz_t multInv(mpz_t a) /* Multiplicative inverse */
{
	//~ A * invA mod N == 1
	mpz_t rop;
	mpz_init(rop);
	mpz_set(rop, euclidian(a));
	return rop;
}

struct point dbl(struct point *p) /* Double operation */
{
	struct point r;
	mpz_inits(r.x, r.y);
	
	if (mpz_cmp_ui(p->y, 0) == 0) // Multiplicative group properties
	{
		mpz_set_ui(r.x,0);
		mpz_set_ui(r.y,0);
	}
	else
	{
		//~ *** Slope ***
		mpz_t s, s0, s1, s2, s3, s4, s6, s6, s7;
		mpz_inits(s, s0, s1, s2, s3, s4, s6, s6, s7);

		//~ s3 = multInv(2*y)
		mpz_mul_ui(s0, p->y, 2);
		mpz_mod(s1, s0, ec.n);
		mpz_set(s2, multInv(s1));
		mpz_mod(s3, s2);

		//~ s6 = 3x² + A
		mpz_powm_ui(s4, p->x, 2, ec.n);
		mpz_mul_ui(s5, s4, 3);
		mpz_add(s6, s5, ec.a);
		
		//~ s = mod(s6 * s3)
		mpz_mul(s7, s6, s3);
		mpz_mod(s, s7, ec.n);

		mpz_clears(s, s0, s1, s2, s3, s4, s6, s6, s7);

		//~ *** Rx ***
		mpz_t rx0, rx1, rx2, rx3, rx4, rx5;
		mpz_inits(rx0, rx1, rx2, rx3, rx4, rx5);

		//~ r.x = s² - (2 * p.x));
		mpz_mul_ui(r0, p->x, 2);
		mpz_mod(r2, r1, ec.n);				
		mpz_powm_ui(r3, s, 2, ec.n);
		mpz_sub(r4, r3,r2);
		mpz_mod(r5, r4, ec.n);
		mpz_set(r.x, r5);

		//~ *** Ry ***
		mpz_t ry0, ry1, ry2, ry3, ry4, ry5;
		mpz_inits(ry0, ry1, ry2, ry3, ry4, ry5);

		//~ r.y = mod(mod(s * mod(p.x - r.x)) - p.y);
		mpz_sub(ry0,p->x, r.x);
		mpz_mod(ry1, ry0, ec.n);
		mpz_mul(ry2, s, ry1);
		mpz_mod(ry3, ry2);
		mpz_sub(ry4, ry3, p->y);
		mpz_mod(ry5, ry4, ec.n);
		mpz_set(r.y, ry5);
	}
	return r;
}

struct point add(struct point *p, struct point *q) /* Add operation */
{
	//~ R = P + Q
	struct coord r;

	//~ *** Slope ***
	mpz_t s, s0, s1, s2, s3, s4, s6, s6
	mpz_inits(s, s0, s1, s2, s3, s4, s6, s6);

	//~ int s = mod(mod(q.y - p.y) * mod(multInv(mod(q.x - p.x))));	// Slope
	mpz_sub(s0, q->y, p->y);
	mpz_mod(s1, s0, ec.n);
	mpz_sub(s2, q->x, p->x);
	mpz_mod(s3, s2, ec.n);
	mpz_mod(s4, multInv(s3), ec.n);
	mpz_mul(s5, s1, s4);
	mpz_mod(s6, s5, ec.n);
	mpz_set(s, s6);

	//~ *** Rx ***
	mpz_t rx0, rx1, rx2, rx3, rx4;
	mpz_inits(rx0, rx1, rx2, rx3, rx4);
	
	//~ r.x = mod(mod(mod(s * s) - p.x) - q.x);
	mpz_powm_ui(rx0, s, 2, ec.n);
	mpz_sub(rx1, rx0, p->x);
	mpz_mod(rx2, rx1, ec.n);
	mpz_sub(rx3, rx2, q->x);
	mpz_mod(rx4, rx3, ec.n);
	mpz_set(r.x, rx4);
	mpz_clears(rx0, rx1, rx2, rx3, rx4);
	
	//~ *** Ry ***
	mpz_t ry0, ry1, ry2, ry3, ry4, ry5;
	mpz_inits(ry0, ry1, ry2, ry3, ry4, ry5);

	//~ r.y = mod(mod(s * mod(p.x - r.x)) - p.y);

	return r;
}

struct point mult( unsigned int k, struct point p) /* Multiplication */
{
	//~ R = kP. Left-to-right binary algorithm.
	struct coord r0 = p;
	struct coord r1 = p;
	int n = k;
	for (int i = n - 2; i > 0; i--)
	{
		if (i > 1)
			r0 = eccDbl(r0);
		else if (i == 1)
			r0 = eccAdd(r1, r0);
	}
	return r0;
}

struct point sub(struct point p, struct point q) /* subtraction */
{
	//~ R = P - Q
	addInv(&q);
	return add(p, q); 
}
