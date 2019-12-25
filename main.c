#include "main.h"

void showPoints()
{
	printf("\nPoints:\n");
	for (int i = 0; i < nPoints; i++) {
		if (isValidPoint(points[i].x, points[i].y)) {
			gmp_printf("%d: (%Zd,%Zd)\n", i, points[i].x, points[i].y);
		} else {
			printf("point %d is not valid\n", i);
		}
	}
	printf("\n");
}

void showPoint(struct coord p)
{
	printf("\nPoint:\n");
	gmp_printf("(%Zd,%Zd) ", p.x, p.y);
	printf("\n");
}

void modn(mpz_t rop, mpz_t a)
{
	mpz_t tmp;
	mpz_init_set(tmp, a);
	if(mpz_cmp_ui(a, 0) < 0) {
		mpz_add(tmp, a, ec.n);
	}
	mpz_mod(rop, tmp, ec.n); // return
	mpz_clear(tmp);
}

void modp(mpz_t rop, mpz_t a)
{
	mpz_t tmp;
	mpz_init_set(tmp, a);
	if(mpz_cmp_ui(a, 0) < 0) {
		mpz_add(tmp, a, ec.n);
	}
	mpz_mod(rop, tmp, ec.n); // return
	mpz_clear(tmp);
}

void eccDbl(struct coord *rop, struct coord p) /* Double operation */
{
	mpz_t s0, s, s1, s2;
	mpz_init(s0);
	// s0 = -p.y
	mpz_neg(s0, p.y);
	modp(s0, s0);
	// if (p.y == -p.y || p.y == 0)
	if (mpz_cmp(p.y, s0) == 0 || mpz_cmp_ui(p.y, 0) == 0) {
		gmp_printf("ret inf, y=0 ou igual -y: (%Zd, %Zd)\n (%Zd, %zd)\n", p.x, p.y);
		mpz_set_ui(rop->x, 0);
		mpz_set_ui(rop->y, 0);
		rop->inf = INF;
	} else {
		// s = mod(mod((p.x * p.x) * 3 + A) * mod(multInv(mod(p.y * 2))))
		mpz_init(s);
		mpz_init(s1);
		mpz_init(s2);
		mpz_mul_ui(s0, p.y, 2); // s0 = p.y * 2
		//gmp_printf("p.y * 2: %Zd\n", s0);
		modp(s1, s0); // s1 = mod(p.y * 2)
		//gmp_printf("mod(p.y * 2): %Zd\n", s1);
		multInv(s0, s1); // s0 = multinv(mod(p.y * 2))

		modp(s1, s0); // s1 = mod(multinv(mod(p.y * 2)))

		mpz_mul(s0, p.x, p.x); // s0 = p.x * p.x
		modp(s0, s0);
		mpz_mul_ui(s2, s0, 3); // s2 = (p.x * p.x) * 3
		modp(s2, s2);
		mpz_add(s0, s2, ec.a); // s0 = (p.x * p.x) * 3 + A
		modp(s2, s0); // s2 = mod((p.x * p.x) * 3 + A)
		mpz_mul(s0, s2, s1); // s0 = mod((p.x * p.x) * 3 + A) * mod(multInv(mod(p.y * 2)))
		// s = mod(mod((p.x * p.x) * 3 + A) * mod(multInv(mod(p.y * 2))))
		modp(s, s0);
		//gmp_printf("slope: %Zd\n", s);

		//~ *** Rx ***
		// rop.x = mod(mod(s * s) - mod(p.x * 2))
		mpz_mul(s0, s, s); // s0 = (s * s)		
		modp(s1, s0); // s1 = mod(s * s)
		mpz_mul_ui(s0, p.x, 2); // s0 = p.x * 2
		modp(s2, s0); // s2 = mod(p.x * 2)
		mpz_sub(s0, s1, s2); // s0 = mod(s * s) - mod(p.x * 2)
		// rop.x = mod(mod(s * s) - mod(p.x * 2)) - return
		modp(rop->x, s0);
		//gmp_printf("novox: %Zd\n", rop->x);
		//~ *** Ry ***
		// rop.y = mod(mod(s * mod(p.x - rop.x)) - p.y)
		mpz_sub(s0, p.x, rop->x); // s0 = p.x - rop.x
		modp(s1, s0); // s1 = mod(p.x - rop.x)
		mpz_mul(s0, s, s1); // s0 = s * mod(p.x - rop.x)
		modp(s1, s0); // s1 = mod(s * mod(p.x - rop.x))
		mpz_sub(s0, s1, p.y); // s0 = mod(s * mod(p.x - rop.x)) - p.y
		// rop.y = mod(mod(s * mod(p.x - rop.x)) - p.y) - return
		modp(rop->y, s0);
		//gmp_printf("novoy: %Zd\n", rop->y);
		rop->inf = 0;

		mpz_clear(s);
		mpz_clear(s1);
		mpz_clear(s2);
	}
	mpz_clear(s0);
}

void eccAdd(struct coord *rop, struct coord p, struct coord q)
{
	//gmp_printf("adding (%Zd, %Zd)\n and: (%Zd, %Zd)\n", p.x, p.y, q.x, q.y);
	//~ r = p + q
	// same x different lines means line parallel with x, so if P and Q are both
	// on the curve, only point left to hit is at infinity
	if (mpz_cmp(p.x, q.x) == 0) {
		if (mpz_cmp(p.y, q.y) != 0) {
			// return
			//gmp_printf("ret inf, mesmo x, dif y: (%Zd, %Zd)\n (%Zd, %zd)\n", p.x, p.y, q.x, q.y);
			mpz_set_ui(rop->x, 0);
			mpz_set_ui(rop->y, 0);
			rop->inf = INF;
			return;
		} else {
			//printf("doubling then\n");
			// same x and y means p + p = 2p
			eccDbl(rop, p);
			return;
		}
	}
	// O + q = q
	if (p.inf == INF) {
		// return
		mpz_set(rop->x, q.x);
		mpz_set(rop->y, q.y);
		rop->inf = q.inf;
		return;
	// p + O = p
	} else if (q.inf == INF) {
		// return
		mpz_set(rop->x, p.x);
		mpz_set(rop->y, p.y);
		rop->inf = p.inf;
		return;
	}

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
	rop->inf = 0;

	mpz_clear(s);
	mpz_clear(s0);
	mpz_clear(s1);
	mpz_clear(s2);
	mpz_clear(s3);
}

void eccSub(struct coord *rop, struct coord p, struct coord q)
{
	mpz_t s0;
	mpz_init(s0);

	struct coord q1;
	mpz_init_set(q1.x, q.x); // q1.x = q.x
	mpz_init(q1.y);
	q1.inf = q.inf;
	// q1.y = mod(q.y * -1);
	mpz_neg(s0, q.y);
	modp(q1.y, s0);
	eccAdd(rop, p, q1);
	mpz_clear(s0);
	mpz_clear(q1.x);
	mpz_clear(q1.y);
}

void mult(struct coord *rop, mpz_t k, struct coord p)
{
	struct coord r0;
	mpz_init_set(r0.x, p.x);
	mpz_init_set(r0.y, p.y);
	r0.inf = p.inf;
	int size = mpz_sizeinbase(k, 2);
	struct coord aux;
	mpz_init(aux.x);
	mpz_init(aux.y);
	aux.inf = 0;
	//printf("size: %d\n", size);
	// printf("size: %d\n", size);
	// Q <- INF
	mpz_set_ui(rop->x, 0);
	mpz_set_ui(rop->y, 0);
	rop->inf = INF;
	for (mp_bitcnt_t i = 0; i < size; i++) {
		//printf("loop: %ld and bit is: %d\n", i, mpz_tstbit(k, i));
		if (mpz_tstbit(k, i) == 1) {
			mpz_set(aux.x, rop->x);
			mpz_set(aux.y, rop->y);
			aux.inf = rop->inf;
			//printf("adding\n");
			eccAdd(rop, aux, r0);
		}
		//gmp_printf("doubling point p: (%Zd, %Zd)\n", r0.x, r0.y);
		mpz_set(aux.x, r0.x);
		mpz_set(aux.y, r0.y);
		aux.inf = r0.inf;
		eccDbl(&r0, aux);
		//gmp_printf("r0: (%Zd, %Zd)\n", r0.x, r0.y);
	}

	mpz_clear(r0.x);
	mpz_clear(r0.y);
}

void euclidian(mpz_t rop, mpz_t a) /* Extended Euclidian algorithm */
{
	mpz_t u, v, g, u1, v1, g1, q, t1, t1t, t2, t2t, t3, t3t;

	mpz_init_set_ui(u,1);
	mpz_init_set_ui(v,0);
	mpz_init_set(g, a);
	mpz_init_set_ui(u1,0);
	mpz_init_set_ui(v1,1);
	mpz_init_set(g1, ec.n);
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

	// alice key pair
	gmp_randseed(state, seed);
	mpz_urandomm(s0, state, ec.n); // s0 = rand % n
	// PEs[i].k = (rand % p) + 1
	mpz_init(alice.k);
	mpz_add_ui(alice.k, s0, 1); //~ Chave privada (k)
	gmp_printf("chave privada a: %Zd\n", alice.k);
	mpz_init(alice.pu.x);
	mpz_init(alice.pu.y);
	alice.pu.inf = 0;
	mult(&alice.pu, alice.k, ec.G);	//~ Chave Pública (pu)
	gmp_printf("pub a: (%Zd, %Zd)\n", alice.pu.x, alice.pu.y);
	// initialize session key
	mpz_init(alice.k_sess);


	// bob key pair
	mpz_add(seed, seed, seed);
	gmp_randseed(state, seed);
	mpz_urandomm(s0, state, ec.n); // s0 = rand % n
	// PEs[i].k = (rand % p) + 1
	mpz_init(bob.k);
	mpz_add_ui(bob.k, s0, 1); //~ Chave privada (k)
	gmp_printf("chave privada b: %Zd\n", bob.k);
	mpz_init(bob.pu.x);
	mpz_init(bob.pu.y);
	bob.pu.inf = 0;
	mult(&bob.pu, bob.k, ec.G);	//~ Chave Pública (pu)
	gmp_printf("pub b: (%Zd, %Zd)\n", bob.pu.x, bob.pu.y);
	// initialize session keys
	mpz_init(bob.k_sess);

	mpz_clear(s0);
	mpz_clear(s1);
	mpz_clear(seed);
	gmp_randclear(state);
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

void findPoints()
{
	mpz_init_set_str(points[0].x, "6B17D1F2E12C4247F8BCE6E563A440F277037D812DEB33A0F4A13945D898C296", 16);
	mpz_init_set_str(points[0].y, "4FE342E2FE1A7F9B8EE7EB4A7C0F9E162BCE33576B315ECECBB6406837BF51F5", 16);
	points[0].inf = 0;
	mpz_init_set_str(points[1].x, "7CF27B188D034F7E8A52380304B51AC3C08969E277F21B35A60B48FC47669978", 16);
	mpz_init_set_str(points[1].y, "07775510DB8ED040293D9AC69F7430DBBA7DADE63CE982299E04B79D227873D1", 16);
	points[1].inf = 0;
	mpz_init_set_str(points[2].x, "5ECBE4D1A6330A44C8F7EF951D4BF165E6C6B721EFADA985FB41661BC6E7FD6C", 16);
	mpz_init_set_str(points[2].y, "8734640C4998FF7E374B06CE1A64A2ECD82AB036384FB83D9A79B127A27D5032", 16);
	points[2].inf = 0;
	mpz_init_set_str(points[3].x, "E2534A3532D08FBBA02DDE659EE62BD0031FE2DB785596EF509302446B030852", 16);
	mpz_init_set_str(points[3].y, "E0F1575A4C633CC719DFEE5FDA862D764EFC96C3F30EE0055C42C23F184ED8C6", 16);
	points[3].inf = 0;
	mpz_init_set_str(points[4].x, "51590B7A515140D2D784C85608668FDFEF8C82FD1F5BE52421554A0DC3D033ED", 16);
	mpz_init_set_str(points[4].y, "E0C17DA8904A727D8AE1BF36BF8A79260D012F00D4D80888D1D0BB44FDA16DA4", 16);
	points[4].inf = 0;
	mpz_init_set_str(points[5].x, "B01A172A76A4602C92D3242CB897DDE3024C740DEBB215B4C6B0AAE93C2291A9", 16);
	mpz_init_set_str(points[5].y, "E85C10743237DAD56FEC0E2DFBA703791C00F7701C7E16BDFD7C48538FC77FE2", 16);
	points[5].inf = 0;
	mpz_init_set_str(points[6].x, "8E533B6FA0BF7B4625BB30667C01FB607EF9F8B8A80FEF5B300628703187B2A3", 16);
	mpz_init_set_str(points[6].y, "73EB1DBDE03318366D069F83A6F5900053C73633CB041B21C55E1A86C1F400B4", 16);
	points[6].inf = 0;
	mpz_init_set_str(points[7].x, "62D9779DBEE9B0534042742D3AB54CADC1D238980FCE97DBB4DD9DC1DB6FB393", 16);
	mpz_init_set_str(points[7].y, "AD5ACCBD91E9D8244FF15D771167CEE0A2ED51F6BBE76A78DA540A6A0F09957E", 16);
	points[7].inf = 0;
	mpz_init_set_str(points[8].x, "EA68D7B6FEDF0B71878938D51D71F8729E0ACB8C2C6DF8B3D79E8A4B90949EE0", 16);
	mpz_init_set_str(points[8].y, "2A2744C972C9FCE787014A964A8EA0C84D714FEAA4DE823FE85A224A4DD048FA", 16);
	points[8].inf = 0;
	mpz_init_set_str(points[9].x, "CEF66D6B2A3A993E591214D1EA223FB545CA6C471C48306E4C36069404C5723F", 16);
	mpz_init_set_str(points[9].y, "878662A229AAAE906E123CDD9D3B4C10590DED29FE751EEECA34BBAA44AF0773", 16);
	points[9].inf = 0;
	mpz_init_set_str(points[10].x, "3ED113B7883B4C590638379DB0C21CDA16742ED0255048BF433391D374BC21D1", 16);
	mpz_init_set_str(points[10].y, "9099209ACCC4C8A224C843AFA4F4C68A090D04DA5E9889DAE2F8EEFCE82A3740", 16);
	points[10].inf = 0;
	mpz_init_set_str(points[11].x, "741DD5BDA817D95E4626537320E5D55179983028B2F82C99D500C5EE8624E3C4", 16);
	mpz_init_set_str(points[11].y, "0770B46A9C385FDC567383554887B1548EEB912C35BA5CA71995FF22CD4481D3", 16);
	points[11].inf = 0;
	mpz_init_set_str(points[12].x, "177C837AE0AC495A61805DF2D85EE2FC792E284B65EAD58A98E15D9D46072C01", 16);
	mpz_init_set_str(points[12].y, "63BB58CD4EBEA558A24091ADB40F4E7226EE14C3A1FB4DF39C43BBE2EFC7BFD8", 16);
	points[12].inf = 0;
	mpz_init_set_str(points[13].x, "54E77A001C3862B97A76647F4336DF3CF126ACBE7A069C5E5709277324D2920B", 16);
	mpz_init_set_str(points[13].y, "F599F1BB29F4317542121F8C05A2E7C37171EA77735090081BA7C82F60D0B375", 16);
	points[13].inf = 0;
	mpz_init_set_str(points[14].x, "F0454DC6971ABAE7ADFB378999888265AE03AF92DE3A0EF163668C63E59B9D5F", 16);
	mpz_init_set_str(points[14].y, "B5B93EE3592E2D1F4E6594E51F9643E62A3B21CE75B5FA3F47E59CDE0D034F36", 16);
	points[14].inf = 0;
	mpz_init_set_str(points[15].x, "76A94D138A6B41858B821C629836315FCD28392EFF6CA038A5EB4787E1277C6E", 16);
	mpz_init_set_str(points[15].y, "A985FE61341F260E6CB0A1B5E11E87208599A0040FC78BAA0E9DDD724B8C5110", 16);
	points[15].inf = 0;
	mpz_init_set_str(points[16].x, "47776904C0F1CC3A9C0984B66F75301A5FA68678F0D64AF8BA1ABCE34738A73E", 16);
	mpz_init_set_str(points[16].y, "AA005EE6B5B957286231856577648E8381B2804428D5733F32F787FF71F1FCDC", 16);
	points[16].inf = 0;
	mpz_init_set_str(points[17].x, "1057E0AB5780F470DEFC9378D1C7C87437BB4C6F9EA55C63D936266DBD781FDA", 16);
	mpz_init_set_str(points[17].y, "F6F1645A15CBE5DC9FA9B7DFD96EE5A7DCC11B5C5EF4F1F78D83B3393C6A45A2", 16);
	points[17].inf = 0;
	mpz_init_set_str(points[18].x, "CB6D2861102C0C25CE39B7C17108C507782C452257884895C1FC7B74AB03ED83", 16);
	mpz_init_set_str(points[18].y, "58D7614B24D9EF515C35E7100D6D6CE4A496716E30FA3E03E39150752BCECDAA", 16);
	points[18].inf = 0;
	mpz_init_set_str(points[19].x, "83A01A9378395BAB9BCD6A0AD03CC56D56E6B19250465A94A234DC4C6B28DA9A", 16);
	mpz_init_set_str(points[19].y, "76E49B6DE2F73234AE6A5EB9D612B75C9F2202BB6923F54FF8240AAA86F640B8", 16);
	points[19].inf = 0;
}

void eccCipher(struct coord *c2, struct coord *c1, int size)
{
	// struct coord aux1, aux2;
	// mpz_init(aux1.x);
	// mpz_init(aux1.y);
	// mpz_init(aux2.x);
	// mpz_init(aux2.y);
	// aux1.inf = 0;
	// aux2.inf = 0;
	int groupAmnt;
	struct coord aux;
	mpz_init(aux.x);
	mpz_init(aux.y);
	// random
	mpz_t seed;
	gmp_randstate_t state;
	gmp_randinit_mt(state);
	mpz_init_set_str(seed, "C49D360886E704936A6678E1139D26B7819F7E90", 16);
	gmp_randseed(state, seed);
	mpz_init(alice.k_sess);
	mpz_urandomm(alice.k_sess, state, ec.n); // chave privada de sessão
	// chave pública de sessão c1 = K = kG	mpz_init(c1.x);
	mpz_init(c1->y);
	mpz_init(c1->x);
	c1->inf = 0;
	mult(c1, alice.k_sess, ec.G);
	groupAmnt = (size % 15 == 0 ? size/15 : size/15 + 1);
	printf("pontos enviados: \n");
	for (int i = 0; i < groupAmnt; i++) {
		mpz_init(c2[i].x);
		mpz_init(c2[i].y);
		c2[i].inf = 0;
		// ponto cifrado 2 é a soma do ponto claro com a multiplicação da chave
		// privada de sessão com a chave pública do destinatário
		mult(&aux, alice.k_sess, bob.pu); // trocar k por chave privada de sessão
		//gmp_printf("cipher: adding (%Zd, %Zd)\n (%Zd, %Zd)\n", points[i % nPoints].x, points[i % nPoints].y, aux.x, aux.y);

		eccAdd(&c2[i], points[i % nPoints], aux);
		gmp_printf("(%Zd, %Zd)\n", points[i % nPoints].x, points[i % nPoints].y);
	}

	printf("Generated ciphered message in points: \n");
	for (int i = 0; i < groupAmnt; i++) {
		gmp_printf("(%Zd,%Zd)\n", c2[i].x, c2[i].y);
	}
	printf("\n");
	mpz_clear(aux.x);
	mpz_clear(aux.y);
	mpz_clear(seed);
	gmp_randclear(state);
}

void eccDecipher(struct coord *d, struct coord *c2, struct coord c1, int size)
{
	int groupAmnt;
	struct coord aux;
	mpz_init(aux.x);
	mpz_init(aux.y);
	groupAmnt = (size % 15 == 0 ? size/15 : size/15 + 1);
	for (int i = 0; i < groupAmnt; i++) {
		mpz_init(d[i].x);
		mpz_init(d[i].y);
		d[i].inf = 0;
		// ponto decifrado é igual a subtração entre o ponto cifrado e o resultado da
		// multiplicação entre a chave privada do destinatário com a chave pública de sessão
		mult(&aux, bob.k, c1);
		eccSub(&d[i], c2[i], aux);
	}

	printf("Generated deciphered message in points: \n");
	for (int i = 0; i < groupAmnt; i++) {
		gmp_printf("(%Zd,%Zd)\n", d[i].x, d[i].y);
	}
	printf("\n");
	mpz_clear(aux.x);
	mpz_clear(aux.y);
}

void main()
{
	mpz_init_set_str(ec.p, "FFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFF", 16);
	mpz_init_set_str(ec.a, "FFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFC", 16);
	mpz_init_set_str(ec.b, "5AC635D8AA3A93E7B3EBBD55769886BC651D06B0CC53B0F63BCE3C3E27D2604B", 16);
	mpz_init_set_str(ec.G.x, "6B17D1F2E12C4247F8BCE6E563A440F277037D812DEB33A0F4A13945D898C296", 16);
	mpz_init_set_str(ec.G.y, "4FE342E2FE1A7F9B8EE7EB4A7C0F9E162BCE33576B315ECECBB6406837BF51F5", 16);
	ec.G.inf = 0;
	mpz_init_set_str(ec.n, "FFFFFFFF00000000FFFFFFFFFFFFFFFFBCE6FAADA7179E84F3B9CAC2FC632551", 16);

	// mpz_init_set_ui(ec.p, 23);
	// mpz_init_set_ui(ec.a, 1);
	// mpz_init_set_ui(ec.b, 4);
	// mpz_init_set_ui(ec.G.x, 1);
	// mpz_init_set_ui(ec.G.y, 12);
	// ec.G.inf = 0;
	// mpz_init_set_ui(ec.n, 3);

	// gmp_printf("ec.p: %Zd\n", ec.p);
	// gmp_printf("ec.a: %Zd\n", ec.a);
	// gmp_printf("ec.b: %Zd\n", ec.b);
	
	findPoints();
	// mpz_t dois;
	// mpz_init_set_ui(dois, 2);
	// struct coord p, q, r1;
	// mpz_init_set_str(p.x, "48439561293906451759052585252797914202762949526041747995844080717082404635286", 10);
	// mpz_init_set_str(p.y, "36134250956749795798585127919587881956611106672985015071877198253568414405109", 10);
	// p.inf=0;
	// mpz_init_set_str(q.x, "48439561293906451759052585252797914202762949526041747995844080717082404635286", 10);
	// mpz_init_set_str(q.y, "36134250956749795798585127919587881956611106672985015071877198253568414405109", 10);
	// q.inf=0;
	// mpz_init(r1.x);
	// mpz_init(r1.y);
	// r1.inf=0;
	// mult(&r1, dois, p);
	// gmp_printf("mult 2 * p: (%Zd, %Zd)\n", r1.x, r1.y);
	// eccDbl(&r1, p);
	// gmp_printf("double p: (%Zd, %Zd)\n", r1.x, r1.y);


	showPoints();
	generateKeys();
	gmp_printf("G:[%Zd,%Zd]\n", ec.G.x, ec.G.y);
	gmp_printf("PU[Alice]:[%Zd,%Zd]\n", alice.pu.x, alice.pu.y);
	gmp_printf("PU[Bob]:[%Zd,%Zd]\n", bob.pu.x, bob.pu.y);

	int messageSize = 32;
	int groupAmnt = (messageSize % 15 == 0 ? messageSize/15 : messageSize/15 + 1);
	// Cifra
	struct coord ciphered[groupAmnt];
	struct coord pub_sess;
	eccCipher(ciphered, &pub_sess, messageSize);

	// Decifra
	struct coord deciphered[groupAmnt];
	eccDecipher(deciphered, ciphered, pub_sess, messageSize);

	mpz_clear(ec.p);
	mpz_clear(ec.a);
	mpz_clear(ec.b);
	mpz_clear(ec.G.x);
	mpz_clear(ec.G.y);
	mpz_clear(ec.n);
}