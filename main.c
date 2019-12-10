#include "main.h"

void showPoints()
{
	printf("\nPoints:\n");
	for (int i = 0; i < pamnt; i++) {
		gmp_printf("%d: (%Zd,%Zd)\n", i, points[i].x, points[i].y);
	}
	printf("\n");
}

void showPoint(struct coord p)
{
	printf("\nPoint:\n");
	gmp_printf("(%Zd,%Zd) ", p.x, p.y);
	printf("\n");
}

void modn(mpz_t res, mpz_t a)
{
	mpz_t rop, tmp;
	mpz_init(rop);	
	mpz_init_set(tmp, a);
	if(mpz_cmp_ui(a, 0) < 0)
		mpz_add(tmp, a, ec.n);
	mpz_mod(rop, tmp, ec.n);
	mpz_set(res, rop);
	mpz_clear(tmp);
	mpz_clear(rop);
}

void modp(mpz_t res, mpz_t a)
{
	mpz_t rop, tmp;
	mpz_init(rop);	
	mpz_init_set(tmp, a);
	if(mpz_cmp_ui(a, 0) < 0)
		mpz_add(tmp, a, ec.p);
	mpz_mod(rop, tmp, ec.p);
	mpz_set(res, rop);
	mpz_clear(tmp);
	mpz_clear(rop);
}

void eccAdd(struct coord *rop, struct coord p, struct coord q)
{
	//~ r = p + q
	//~ *** Slope ***
	mpz_t s, s0, s1, s2, s3;
	mpz_init(s);
	mpz_init(s0);
	mpz_init(s1);
	mpz_init(s2);
	mpz_init(s3);

	//gmp_printf("sum of p: %Zd, %Zd\nand q: %Zd, %Zd\n", p.x, p.y, q.x, q.y);

	// s = mod(mod(q.y - p.y) * mod(multInv(mod(q.x - p.x))));	// Slope
	// s0 = q.y - p.y
	mpz_sub(s0, q.y, p.y);
	// s1 = mod(q.y - p.y)
	modp(s1, s0);
	//gmp_printf("mod(q.y - p.y): %Zd\n", s1);
	// s0 = q.x - p.x
	mpz_sub(s0, q.x, p.x);
	// s2 = mod(q.x - p.x)
	modp(s2, s0);
	// s0 = multInv(mod(q.x - p.x))
	multInv(s0, s2);
	// s3 = mod(multInv(mod(q.x - p.x)))
	modp(s3, s0);
	// s2 = mod(1.y - p.y) * mod(multInv(mod(q.x - p.x)))
	mpz_mul(s2, s1, s3);
	// s = mod(mod(1.y - p.y) * mod(multInv(mod(q.x - p.x))))
	modp(s, s2);

	//~ *** Rx ***
	// rop.x = mod(mod(s * s) - p.x) - q.x)
	// s0 = s * s
	mpz_mul(s0, s, s);
	// s1 = mod(s * s)
	modp(s1, s0);
	// s0 = mod(s * s) - p.x
	mpz_sub(s0, s1, p.x);
	// s1 = mod(mod(s * s) - p.x)
	modp(s1, s0);
	// s0 = mod(mod(s * s) - p.x) - q.x;
	mpz_sub(s0, s1, q.x);
	// rop.x = mod(mod(s * s) - p.x) - q.x)
	modp(rop->x, s0);
	
	//~ *** Ry ***
	// rop.y = mod(mod(s * mod(p.x - rop.x)) - p.y)
	// rop.y = mod(s0)
	// s0 = p.x - rop.x
	mpz_sub(s0, p.x, rop->x);
	// s1 = mod(p.x - rop.x)
	modp(s1, s0);
	// s0 = s * mod(p.x - rop.x)
	mpz_mul(s0, s, s1);
	// s1 = mod(s * mod(p.x - rop.x))
	modp(s1, s0);
	// s0 = mod(s * mod(p.x - rop.x)) - p.y
	mpz_sub(s0, s1, p.y);
	// rop.y = mod(mod(s * mod(p.x - rop.x)) - p.y)
	modp(rop->y, s0);

	// return
	modp(rop->x, rop->x);
	modp(rop->y, rop->y);

	mpz_clear(s);
	mpz_clear(s0);
	mpz_clear(s1);
	mpz_clear(s2);
	mpz_clear(s3);
}

void eccDbl(struct coord *rop, struct coord p) /* Double operation */
{
	mpz_t s0, s, s1, s2;
	mpz_init(s0);
	// s0 = -p.y
	mpz_neg(s0, p.y);
	// if (p.y == -p.y || p.y == 0)
	if (mpz_cmp(p.y, s0) == 0 || mpz_cmp_ui(p.y, 0) == 0) {
		mpz_set_ui(rop->x, 0);
		mpz_set_ui(rop->y, 0);
	} else {
		//~ *** Slope ***
		// s = mod(mod((p.x * p.x) * 3 + A) * mod(multInv(mod(p.y * 2))))
		mpz_init(s);
		mpz_init(s1);
		mpz_init(s2);
		mpz_mul_ui(s0, p.y, 2); // s0 = p.y * 2
		modp(s1, s0); // s1 = mod(p.y * 2)
		multInv(s0, s1); // s0 = multinv(mod(p.y * 2))
		modp(s1, s0); // s1 = mod(multinv(mod(p.y * 2)))
		mpz_mul(s0, p.x, p.x); // s0 = p.x * p.x
		mpz_mul_ui(s2, s0, 3); // s2 = (p.x * p.x) * 3
		mpz_add(s0, s2, ec.a); // s0 = (p.x * p.x) * 3 + A
		modp(s2, s0); // s2 = mod((p.x * p.x) * 3 + A)
		mpz_mul(s0, s2, s1); // s0 = mod((p.x * p.x) * 3 + A) * mod(multInv(mod(p.y * 2)))
		// s = mod(mod((p.x * p.x) * 3 + A) * mod(multInv(mod(p.y * 2))))
		modp(s, s0);

		//~ *** Rx ***
		// rop.x = mod(mod(s * s) - mod(p.x * 2))
		mpz_mul(s0, s, s); // s0 = (s * s)		
		modp(s1, s0); // s1 = mod(s * s)
		mpz_mul_ui(s0, p.x, 2); // s0 = p.x * 2
		modp(s2, s0); // s2 = mod(p.x * 2)
		mpz_sub(s0, s1, s2); // s0 = mod(s * s) - mod(p.x * 2)
		// rop.x = mod(mod(s * s) - mod(p.x * 2)) - return
		modp(rop->x, s0);

		//~ *** Ry ***
		// rop.y = mod(mod(s * mod(p.x - rop.x)) - p.y)
		mpz_sub(s0, p.x, rop->x); // s0 = p.x - rop.x
		modp(s1, s0); // s1 = mod(p.x - rop.x)
		mpz_mul(s0, s, s1); // s0 = s * mod(p.x - rop.x)
		modp(s1, s0); // s1 = mod(s * mod(p.x - rop.x))
		mpz_sub(s0, s1, p.y); // s0 = mod(s * mod(p.x - rop.x)) - p.y
		// rop.y = mod(mod(s * mod(p.x - rop.x)) - p.y) - return
		modp(rop->y, s0);

		mpz_clear(s);
		mpz_clear(s1);
		mpz_clear(s2);
	}
	mpz_clear(s0);
}

void eccSub(struct coord *rop, struct coord p, struct coord q)
{
	mpz_t s0;
	mpz_init(s0);

	struct coord q1;
	mpz_init_set(q1.x, q.x); // q1.x = q.x
	mpz_init(q1.y);
	// q1.y = mod(q.y * -1);
	mpz_neg(s0, q.y);
	modp(q1.y, s0);
	eccAdd(rop, p, q1);
	mpz_clear(s0);
	mpz_clear(q1.x);
	mpz_clear(q1.y);
}

//~ Left-to-right binary algorithm
void mult(struct coord *rop, mpz_t k, struct coord p)
{
	//~ output: q = kP
	struct coord r0;
	mpz_init_set(r0.x, p.x);
	mpz_init_set(r0.y, p.y);
	mpz_t i;
	mpz_init(i);
	// n = k - 2
	for (mpz_sub_ui(i, k, 2); mpz_cmp_ui(i, 0) > 0; mpz_sub_ui(i, i, 1)) {
		eccDbl(rop, *rop); // - return
	}
	eccAdd(rop, r0, *rop); // - return
	printf("\n");
	mpz_clear(r0.x);
	mpz_clear(r0.y);
	mpz_clear(i);
}

void euclidian(mpz_t rop, mpz_t a) /* Extended Euclidian algorithm */
{
	mpz_t u, v, g, u1, v1, g1, q, t1, t1t, t2, t2t, t3, t3t;

	mpz_init_set_ui(u,1);
	mpz_init_set_ui(v,0);
	mpz_init_set(g, a);
	mpz_init_set_ui(u1,0);
	mpz_init_set_ui(v1,1);
	mpz_init_set(g1, ec.p);
	mpz_init(q);
	mpz_init(t1);
	mpz_init(t1t);
	mpz_init(t2);
	mpz_init(t2t);
	mpz_init(t3);
	mpz_init(t3t);

	//~ while (g1 != 0)
	while(mpz_cmp_ui(g1, 0) != 0)
	{
		//~ q = g / g1;
		mpz_fdiv_q(q, g, g1); 
		
		//~ t1 = u - q * u1;
		mpz_mul(t1t, q, u1);
		mpz_sub(t1, u, t1t);
		
		//~ t2 = v - q * v1;
		mpz_mul(t2t, q, v1);
		mpz_sub(t2, v, t2t);

		//~ t3 = g - q * g1;
		mpz_mul(t3t, q, g1);
		mpz_sub(t3, g, t3t);

		mpz_set(u, u1);
		mpz_set(v, v1);
		mpz_set(g, g1);
		mpz_set(u1, t1);
		mpz_set(v1, t2);
		mpz_set(g1, t3);
	}
	
	//return
	modp(rop, u);

	mpz_clear(u);
	mpz_clear(v);
	mpz_clear(g);
	mpz_clear(u1);
	mpz_clear(v1);
	mpz_clear(g1);
	mpz_clear(q);
	mpz_clear(t1);
	mpz_clear(t1t);
	mpz_clear(t2);
	mpz_clear(t2t);
	mpz_clear(t3);
	mpz_clear(t3t);
}


void multInv(mpz_t rop, mpz_t a) /* Multiplicative inverse */
{
	//~ A * invA mod P == 1
	// return
	euclidian(rop, a);
}

//~ Gerar as chaves públicas e privadas
void generateKeys()
{
	mpz_t s0, s1, seed;
	gmp_randstate_t state;
	gmp_randinit_mt(state);
	mpz_init_set_str(seed, "C49D360886E704936A6678E1139D26B7819F7E90", 16);
	mpz_init(s0);
	mpz_init(s1);
	for (int i = 0; i < N; i++)
	{
		gmp_randseed(state, seed);
		mpz_urandomm(s0, state, ec.p); // s0 = rand % p
		// PEs[i].k = (rand % p) + 1
		mpz_init(PEs[i].k);
		mpz_add_ui(PEs[i].k, s0, 1); //~ Chave privada (k)
		gmp_printf("chave privada %d: %Zd\n", i, PEs[i].k);
		mpz_init(PEs[i].pu.x);
		mpz_init(PEs[i].pu.y);
		mult(&PEs[i].pu, PEs[i].k, ec.G);	//~ Chave Pública (pu)
		gmp_printf("pub %d: %Zd, %Zd\n", i, PEs[i].pu.x, PEs[i].pu.y);
		mpz_add(seed, seed, seed);
	}
	mpz_clear(s0);
	mpz_clear(s1);
	mpz_clear(seed);
	gmp_randclear(state);
}

void generateMessage(int size, struct coord *m, char *plaintext)
{
	for (int i = 0; i < size; i++)
	{
		int index = (rand() % 25);
		plaintext[i] = alpha[index];
		mpz_init_set(m[i].x, points[index].x);
		mpz_init_set(m[i].y, points[index].y);
	}
	printf("Generated message in chat: %s\n", plaintext);
	printf("Generated message in points: ");
	for (int i = 0; i < size; i++) {
		gmp_printf("(%Zd,%Zd) ", m[i].x, m[i].y);
	}
	printf("\n");
}

void eccCipher(struct pe a, struct pe b, struct coord *m, struct coord *c, int s)
{
	struct coord aux;
	mpz_init(aux.x);
	mpz_init(aux.y);
	for (int i = 0; i < s; i++) {
		mpz_init(c[i].x);
		mpz_init(c[i].y);
		// ponto cifrado é igual a soma do ponto claro e o resultado da multiplicação
		// entre a chave privada do remetente com a chave pública do destinatário
		mult(&aux, a.k, b.pu);
		eccAdd(&c[i], m[i], aux);
	}

	printf("Generated ciphered message in points: \n");
	for (int i = 0; i < s; i++) {
		gmp_printf("(%Zd,%Zd) ", c[i].x, c[i].y);
	}
	printf("\n");
	mpz_clear(aux.x);
	mpz_clear(aux.y);
}

void eccDecipher(struct pe a, struct pe b, struct coord *d, struct coord *c, int s)
{
	struct coord aux;
	mpz_init(aux.x);
	mpz_init(aux.y);
	for (int i = 0; i < s; i++) {
		mpz_init(d[i].x);
		mpz_init(d[i].y);
		// ponto decifrado é igual a subtração do ponto cifrado e o resultado da multiplicação
		// entre a chave privada do destinatário com a chave pública do remetente
		mult(&aux, b.k, a.pu);
		eccSub(&d[i], c[i], aux);
	}

	printf("Generated deciphered message in points: \n");
	for (int i = 0; i < s; i++) {
		gmp_printf("(%Zd,%Zd) ", d[i].x, d[i].y);
	}
	printf("\n");
	mpz_clear(aux.x);
	mpz_clear(aux.y);
}

//~ Verifica a validade do ponto
bool isValidPoint(mpz_t x, mpz_t y)
{
	mpz_t op1, op2, s0, s1, s2;
	mpz_init(op1);
	mpz_init(op2);
	mpz_init(s0);
	mpz_init(s1);
	mpz_init(s2);
	// gmp_printf("checking: (x: %Zd, y: %Zd)\n", x, y);
	// op1 = mod(y * y)
	mpz_mul(s0, y, y);// s0 = y * y
	modp(op1, s0); // op1 = mod(y * y)
	// gmp_printf("op1: %Zd\n", op1);
	// op2 = mod(mod(x * x * x) + mod(A * x) + B)
	mpz_mul(s0, x, x); // s0 = x * x
	mpz_mul(s1, x, s0); // s1 = x * x * x
	modp(s0, s1); // s0 = mod(x * x * x)
	// gmp_printf("mod(x * x * x): %Zd\n", s0);
	mpz_mul(s1, ec.a, x); // s1 = A * x
	modp(s2, s1); // s2 = mod(A * x)
	// gmp_printf("mod(A * x): %Zd\n", s2);
	mpz_add(s1, s2, ec.b); // s1 = mod(A * x) + B
	mpz_add(s2, s0, s1); // s2 = mod(x * x * x) + mod(A * x) + B
	modp(op2, s2); // op2 = mod(mod(x * x * x) + mod(A * x) + B)
	// gmp_printf("op2: %Zd\n", op2);
	if (mpz_cmp(op1, op2) == 0) {
		mpz_clear(op1);
		mpz_clear(op2);
		mpz_clear(s0);
		mpz_clear(s1);
		mpz_clear(s2);
		return true;
	}
	mpz_clear(op1);
	mpz_clear(op2);
	mpz_clear(s0);
	mpz_clear(s1);
	mpz_clear(s2);
	return false;
}

bool isnotduplicate(struct coord c, int k)
{
	for (int i = 0; i < k; i++) {
		if (mpz_cmp(points[i].x, c.x) == 0 && mpz_cmp(points[i].y, c.y) == 0) {
			return false;
		}
	}
	return true;
}

void findPoints()
{
	struct coord test1, test2, res1, res2, res3, res4;
	mpz_init(test1.x);
	mpz_init(test1.y);
	mpz_init(test2.x);
	mpz_init(test2.y);
	mpz_init(res1.x);
	mpz_init(res1.y);
	mpz_init(res2.x);
	mpz_init(res2.y);
	mpz_init(res3.x);
	mpz_init(res3.y);
	mpz_init(res4.x);
	mpz_init(res4.y);
	int k = 2;
	mpz_t dois;
	mpz_init_set_ui(dois, 2);
	mpz_init_set_str(points[0].x, "6B17D1F2E12C4247F8BCE6E563A440F277037D812DEB33A0F4A13945D898C296", 16);
	mpz_init_set_str(points[0].y, "4FE342E2FE1A7F9B8EE7EB4A7C0F9E162BCE33576B315ECECBB6406837BF51F5", 16);
	mpz_init_set_str(points[1].x, "56515219790691171413109057904011688695424810155802929973526481321309856242040", 10);
	mpz_init_set_str(points[1].y, "03377031843712258259223711451491452598088675519751548567112458094635497583569", 10);
	for (int i = 2; k < pamnt && i < pamnt; i++) {
		mpz_set(test1.x, points[i - 2].x);
		mpz_set(test1.y, points[i - 2].y);
		mpz_set(test2.x, points[i - 1].x);
		mpz_set(test2.y, points[i - 1].y);
		eccAdd(&res1, points[i - 2], points[i - 1]);
		eccSub(&res2, points[i - 2], points[i - 1]);
		eccDbl(&res3, points[i - 1]);
		mult(&res4, dois, points[i - 1]);
		if (isValidPoint(res1.x, res1.y) && isnotduplicate(res1, k)) {
			mpz_init_set(points[k].x, res1.x);
			mpz_init_set(points[k].y, res1.y);
			k++;
			gmp_printf("point %d found! x, y: %Zd, %Zd\n", k, points[k-1].x, points[k-1].x);
		}
		if (isValidPoint(res2.x, res2.y) && isnotduplicate(res1, k)) {
			mpz_init_set(points[k].x, res2.x);
			mpz_init_set(points[k].y, res2.y);
			k++;
			gmp_printf("point %d found! x, y: %Zd, %Zd\n", k, points[k-1].x, points[k-1].x);
		}
		if (isValidPoint(res3.x, res3.y) && isnotduplicate(res1, k)) {
			mpz_init_set(points[k].x, res3.x);
			mpz_init_set(points[k].y, res3.y);
			k++;
			gmp_printf("point %d found! x, y: %Zd, %Zd\n", k, points[k-1].x, points[k-1].x);
		}
		if (isValidPoint(res4.x, res4.y) && isnotduplicate(res1, k)) {
			mpz_init_set(points[k].x, res4.x);
			mpz_init_set(points[k].y, res4.y);
			k++;
			gmp_printf("point %d found! x, y: %Zd, %Zd\n", k, points[k-1].x, points[k-1].x);
		}
	}
	printf("found %d points!\n", k);
}

//~ Calcular o número de pontos na curva
void numberPoints(mpz_t res)
{
	mpf_t fp, aux1, aux2;;
	mpf_set_default_prec(512);
  	mpf_init_set_str(fp, "115792089210356248762697446949407573530086143415290314195533631308867097853951.0", 10);
	mpf_init(aux1);
	mpf_init(aux2);
	// npoints = floor(p + 1 + 2 * sqrt(p));
	//gmp_printf("fp                        : %Ff\n", fp);
	// aux1 = sqrt(p)
	mpf_sqrt(aux1, fp);
	//gmp_printf("sqrt(p)                   : %Ff\n", aux1);
	// aux2 = 2 * sqrt(p)
	mpf_mul_ui(aux2, aux1, 2);
	//gmp_printf("2 * sqrt(p)               : %Ff\n", aux2);
	// aux1 = 1 + 2 * sqrt(p)
	mpf_add_ui(aux1, aux2, 1);
	//gmp_printf("1 + 2 * sqrt(p)           : %Ff\n", aux1);
	// aux2 = p + 1 + 2 * sqrt(p)
	mpf_add(aux2, fp, aux1);
	//gmp_printf("p + 1 + 2 * sqrt(p)       : %Ff\n", aux2);
	// aux1 = floor(p + 1 + 2 * sqrt(p))
	mpf_floor(aux1, aux2);
	//gmp_printf("floor(p + 1 + 2 * sqrt(p)): %Ff\n", aux1);
	mpz_set_f(res, aux1); // returning result
	//gmp_printf("npoints: %Zd\n", res);
	mpf_clear(aux1);
	mpf_clear(aux2);
	mpf_clear(fp);
}

//~ Cálculo de todos pontos da curva do ivan
void calculatePoints()
{
	mpz_t i, j;
	mpz_init(nPoints);
	numberPoints(nPoints);
	gmp_printf("nPoints: %Zd\n", nPoints);
	int k = 0;
	// fazer o break para testar a função isvalidpoint
	// #pragma omp parallel for shared(k) collapse(2) private(j) num_threads(3)
	for (mpz_init_set_str(i, "03377031843712258259223711451491452598088675519751548567112458094635497583569", 10); mpz_cmp(i, ec.p) < 0 && k < pamnt; mpz_add_ui(i, i, 1)) {
		for (mpz_init_set_str(j, "03377031843712258259223711451491452598088675519751548567112458094635497583569", 10); mpz_cmp(j, ec.p) < 0 && k < pamnt; mpz_add_ui(j, j, 1)) {
			if (isValidPoint(i, j)) {
				mpz_init_set(points[k].x, i);
				mpz_init_set(points[k].y, j);
				k++;
				gmp_printf("found (i: %Zd, j: %Zd)\n", i, j);
			}
		}
	}
	printf("k: %d\n", k);
	// printf("%d generated points + Infinity Point (0,1,0)\n", k);
}

//~ It's show time
void main()
{
	//double x = pow(2,192) - pow(2,32) - pow(2,12) - pow(2,8) - pow(2,7) - pow(2,6) - pow(2,3) - 1;
	//printf("%lf",x);
	mpz_init_set_str(ec.p, "FFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFF", 16);
	mpz_init_set_str(ec.a, "FFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFC", 16);
	mpz_init_set_str(ec.b, "5AC635D8AA3A93E7B3EBBD55769886BC651D06B0CC53B0F63BCE3C3E27D2604B", 16);
	mpz_init_set_str(ec.G.x, "6B17D1F2E12C4247F8BCE6E563A440F277037D812DEB33A0F4A13945D898C296", 16);
	mpz_init_set_str(ec.G.y, "4FE342E2FE1A7F9B8EE7EB4A7C0F9E162BCE33576B315ECECBB6406837BF51F5", 16);
	mpz_init_set_str(ec.n, "FFFFFFFF00000000FFFFFFFFFFFFFFFFBCE6FAADA7179E84F3B9CAC2FC632551", 16);

	// mpz_init_set_ui(ec.p, 23);
	// mpz_init_set_ui(ec.a, 1);
	// mpz_init_set_ui(ec.b, 4);
	// mpz_init_set_ui(ec.G.x, 1);
	// mpz_init_set_ui(ec.G.y, 12);
	// mpz_init_set_ui(ec.n, 3);

	gmp_printf("ec.p: %Zd\n", ec.p);
	gmp_printf("ec.a: %Zd\n", ec.a);
	gmp_printf("ec.b: %Zd\n", ec.b);
	
	findPoints();
	showPoints();
	generateKeys();
	gmp_printf("G:[%Zd,%Zd]\n", ec.G.x, ec.G.y);
	gmp_printf("PU[Alice]:[%Zd,%Zd]\n", PEs[ALICE].pu.x, PEs[ALICE].pu.y);
	gmp_printf("PU[Bob]:[%Zd,%Zd]\n", PEs[BOB].pu.x, PEs[BOB].pu.y);
	
	// struct coord teste;
	// mpz_init(teste.x);
	// mpz_init(teste.y);
	// mpz_init_set_str(teste.x, "56515219790691171413109057904011688695424810155802929973526481321309856242040", 10);
	// mpz_init_set_str(teste.y, "03377031843712258259223711451491452598088675519751548567112458094635497583569", 10);
	// eccDbl(&teste, ec.G);
	// gmp_printf("teste.x: %Zd, teste.y: %Zd\n", teste.x, teste.y);
	// if (isValidPoint(teste.x, teste.y)) {
	// 	printf("deu caralho22!\n");
	// }

	int messageSize = 32;		// bytes
	char *plaintext = malloc(sizeof(char) * messageSize);
	struct coord message[messageSize];

	generateMessage(messageSize, message, plaintext);

	// Cifra
	struct coord ciphered[messageSize];
	eccCipher(PEs[ALICE], PEs[BOB], message, ciphered, messageSize);

	// Decifra
	struct coord deciphered[messageSize];
	eccDecipher(PEs[ALICE], PEs[BOB], deciphered, ciphered, messageSize);

	// // Assina  
	// struct coord signature = sign(PEs[ALICE].k, plaintext, messageSize);

	// // Verifica a assinatura
	// if (checkSignature(PEs[ALICE].pu, signature, plaintext))
	// 	printf("Assinatura válida!");
	// else
	// 	printf("Assinatura inválida!");

	// printf("\n\n");
	// points = NULL;

	// mpz_clear(nPoints);
	// mpz_clear(ec.p);
	// mpz_clear(ec.a);
	// mpz_clear(ec.b);
	// mpz_clear(ec.G.x);
	// mpz_clear(ec.G.y);
	// mpz_clear(ec.n);
}
