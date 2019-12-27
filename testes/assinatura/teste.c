#include "teste.h"

double timestamp(void)
{
    struct timeval tp;
    gettimeofday(&tp, NULL);

    return ((double)(tp.tv_sec * 1000.0 + tp.tv_usec / 1000.0));
}

inline void modn(mpz_t rop, mpz_t a)
{
	mpz_t tmp;
	mpz_init_set(tmp, a);
	if(mpz_cmp_ui(a, 0) < 0) {
		mpz_add(tmp, a, ec.n);
	}
	mpz_mod(rop, tmp, ec.n);
	mpz_clear(tmp);
}

inline void modp(mpz_t rop, mpz_t a)
{
	mpz_t tmp;
	mpz_init_set(tmp, a);
	if(mpz_cmp_ui(a, 0) < 0) {
		mpz_add(tmp, a, ec.p);
	}
	mpz_mod(rop, tmp, ec.p);
	mpz_clear(tmp);
}

inline void eccDbl(struct coord *rop, struct coord p) /* Double operation */
{
	mpz_t s0, s, s1, s2;
	mpz_init(s0);
	// s0 = -p.y
	mpz_neg(s0, p.y);
	modp(s0, s0);
	// if (p.y == -p.y || p.y == 0)
	if (mpz_cmp(p.y, s0) == 0 || mpz_cmp_ui(p.y, 0) == 0) {
		// gmp_printf("ret inf, y=0 ou igual -y: (%Zd, %Zd)\n (%Zd, %zd)\n", p.x, p.y);
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
		multInv(s0, s1, ec.p); // s0 = multinv(mod(p.y * 2))
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

	mpz_t s, s0, s1, s2;
	mpz_init(s);
	mpz_init(s0);
	mpz_init(s1);
	mpz_init(s2);

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
	multInv(s0, s2, ec.p);
	// s0 = mod(multInv(mod(q.x - p.x)))
	modp(s0, s0);
	// s2 = mod(1.y - p.y) * mod(multInv(mod(q.x - p.x)))
	mpz_mul(s2, s1, s0);
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
	// rop.x = mod(mod(s * s) - p.x) - q.x) - return
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
	// rop.y = mod(mod(s * mod(p.x - rop.x)) - p.y) - return
	modp(rop->y, s0);
	rop->inf = 0;

	modp(rop->x, rop->x);
	modp(rop->y, rop->y);

	mpz_clear(s);
	mpz_clear(s0);
	mpz_clear(s1);
	mpz_clear(s2);
}

inline void eccSub(struct coord *rop, struct coord p, struct coord q)
{
	mpz_t s0;
	mpz_init(s0);

	struct coord q1;
	mpz_init_set(q1.x, q.x); // q1.x = q.x
	mpz_init(q1.y);
	// q1.y = mod(q.y * -1);
	mpz_neg(s0, q.y);
	q1.inf = q.inf;
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
	mpz_clear(aux.x);
	mpz_clear(aux.y);
}

inline void euclidian(mpz_t rop, mpz_t a, mpz_t md) /* Extended Euclidian algorithm */	
{	
	mpz_t u, g, u1, g1, q, t1, t3, aux;	
	mpz_init_set_ui(u,1);		
	mpz_init_set(g, a);	
	mpz_init_set_ui(u1,0);	
	mpz_init_set(g1, md);	
	mpz_init(q);	
	mpz_init(t1);		
	mpz_init(t3);	
	mpz_init(aux);
	//~ while (g1 != 0)	
	while(mpz_cmp_ui(g1, 0) != 0)	
	{	
		//~ q = g / g1;	
		mpz_fdiv_q(q, g, g1); 	
			
		//~ t1 = u - q * u1;	
		mpz_mul(aux, q, u1);	
		mpz_sub(t1, u, aux);	
		mpz_set(u, u1);
		mpz_set(u1, t1);
			
		//~ t2 = v - q * v1;	
		// mpz_mul(aux, q, v1);	
		// mpz_sub(t2, v, aux);	
		// mpz_set(v, v1);	
		// mpz_set(v1, t2);	

		//~ t3 = g - q * g1;	
		mpz_mul(aux, q, g1);	
		mpz_sub(t3, g, aux);
		mpz_set(g, g1);	
		mpz_set(g1, t3);	
	}
	
	//return

	mpz_set(rop, u);

	mpz_clear(u);
	mpz_clear(g);
	mpz_clear(u1);
	mpz_clear(g1);
	mpz_clear(q);
	mpz_clear(t1);	
	mpz_clear(t3);
	mpz_clear(aux);
}

inline void multInv(mpz_t rop, mpz_t a, mpz_t md) /* Multiplicative inverse */
{
	//~ A * invA mod P == 1
	// return
	euclidian(rop, a, md);
}

//~ Gerar as chaves públicas e privadas
inline void generateKeys()
{
	mpz_t s0, s1;
	mpz_init(s0);
	mpz_init(s1);

	// alice key pair
	mpz_urandomm(s0, state, ec.p); // s0 = rand % n
	// PEs[i].k = (rand % p) + 1
	mpz_init(alice.k);
	mpz_add_ui(alice.k, s0, 1); //~ Chave privada (k)
	// gmp_printf("chave privada a: %Zd\n", alice.k);
	mpz_init(alice.pu.x);
	mpz_init(alice.pu.y);
	alice.pu.inf = 0;
	mult(&alice.pu, alice.k, ec.G);	//~ Chave Pública (pu)
	// gmp_printf("pub a: (%Zd, %Zd)\n", alice.pu.x, alice.pu.y);
	// initialize session key
	mpz_init(alice.k_sess);

	// bob key pair
	mpz_urandomm(s0, state, ec.p); // s0 = rand % n
	// PEs[i].k = (rand % p) + 1
	mpz_init(bob.k);
	mpz_add_ui(bob.k, s0, 1); //~ Chave privada (k)
	// gmp_printf("chave privada b: %Zd\n", bob.k);
	mpz_init(bob.pu.x);
	mpz_init(bob.pu.y);
	bob.pu.inf = 0;
	mult(&bob.pu, bob.k, ec.G);	//~ Chave Pública (pu)
	// gmp_printf("pub b: (%Zd, %Zd)\n", bob.pu.x, bob.pu.y);
	// initialize session keys
	mpz_init(bob.k_sess);

	mpz_clear(s0);
	mpz_clear(s1);
}

inline void findPoints()
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

void genMessages()
{
	for (int i = 0; i < 20; i++) {
		mpz_urandomm(mess[i], state, ec.p);
		//gmp_printf("%Zd\n", mess[i]);
	}
}

inline void eccCipher(struct coord *c2, struct coord *c1, int size)
{
	// struct coord aux1, aux2;
	// mpz_init(aux1.x);
	// mpz_init(aux1.y);
	// mpz_init(aux2.x);
	// mpz_init(aux2.y);
	// aux1.inf = 0;
	// aux2.inf = 0;
	struct coord aux;
	mpz_init(aux.x);
	mpz_init(aux.y);
	mpz_urandomm(alice.k_sess, state, ec.p); // chave privada de sessão
	// chave pública de sessão c1 = K = kG	mpz_init(c1.x);
	mpz_init(c1->y);
	mpz_init(c1->x);
	c1->inf = 0;
	mult(c1, alice.k_sess, ec.G);
	//printf("pontos enviados: \n");
	for (int i = 0; i < size; i++) {
		mpz_init(c2[i].x);
		mpz_init(c2[i].y);
		c2[i].inf = 0;
		// ponto cifrado 2 é a soma do ponto claro com a multiplicação da chave
		// privada de sessão com a chave pública do destinatário
		mult(&aux, alice.k_sess, bob.pu); // trocar k por chave privada de sessão
		//gmp_printf("cipher: adding (%Zd, %Zd)\n (%Zd, %Zd)\n", points[i % nPoints].x, points[i % nPoints].y, aux.x, aux.y);

		eccAdd(&c2[i], points[i % nPoints], aux);
		// gmp_printf("(%Zd, %Zd)\n", points[i % nPoints].x, points[i % nPoints].y);
	}

	// printf("Generated ciphered message in points: \n");
	// for (int i = 0; i < size; i++) {
	// 	gmp_printf("(%Zd,%Zd) ", c2[i].x, c2[i].y);
	// }
	// printf("\n");
	mpz_clear(aux.x);
	mpz_clear(aux.y);
}

inline void eccDecipher(struct coord *d, struct coord *c2, struct coord c1, int size)
{
	struct coord aux;
	mpz_init(aux.x);
	mpz_init(aux.y);
	for (int i = 0; i < size; i++) {
		mpz_init(d[i].x);
		mpz_init(d[i].y);
		d[i].inf = 0;
		// ponto decifrado é igual a subtração entre o ponto cifrado e o resultado da
		// multiplicação entre a chave privada do destinatário com a chave pública de sessão
		mult(&aux, bob.k, c1);
		eccSub(&d[i], c2[i], aux);
	}

	// printf("Generated deciphered message in points: \n");
	// for (int i = 0; i < size; i++) {
	// 	gmp_printf("(%Zd,%Zd) ", d[i].x, d[i].y);
	// }
	// printf("\n");
	mpz_clear(aux.x);
	mpz_clear(aux.y);
}

void print_hash(unsigned char hash[]) {
    int idx;
    for (idx = 0; idx < 32; idx++)
        printf("%02x", hash[idx]);
    printf("\n");
}

inline void hash(mpz_t rop, mpz_t m)
{
	// hashes the hexadecimal (in lower case) representation of m
	// returns the decimal representation of the hash in rop
    int res;
	int size = mpz_sizeinbase(m, 2);
	// gmp_printf("\n - m:%Zd (%d bytes)\n", m, size);
    unsigned char hash[32];
	unsigned char *message = malloc(sizeof(unsigned char) * size + 1);
	char txt[64];
	// transform mpz_t m in string message
	FILE *temp = fopen("num.temp", "w+");
	mpz_out_str(temp, 16, m);
	fseek(temp, 0L, SEEK_SET);
	fgets (message, 65, temp);
	fclose(temp);

    SHA256_CTX ctx;
	sha256_init(&ctx);
	sha256_update(&ctx, message, strlen(message));
	sha256_final(&ctx, hash, &res, size, 23);
	//print_hash(hash);
	// convert unsigned char string to char string
	int a = 0;
	for (int i = 0; i < 32; i++) {
		sprintf(&txt[a], "%02x", hash[i]);
		a+=2;
	}
	//printf("txt: %s\n", txt);
    // write hash in rop
	mpz_set_str(rop, txt, 16);
	//gmp_printf("rop: %Zd\n", rop);
	modp(rop, rop);
	//gmp_printf("mod(rop): %Zd\n", rop);
}

inline void generateCertificates()
{
	// Chaves e IDs
	mpz_init_set_str(alice.id.x, "25992258928441372838190306754026437895674821091948799811836327659484854798868", 10);
	mpz_init_set_str(alice.id.y, "39227732877500211262332505102593161927141165977758306466170574199348208861497", 10);
	alice.id.inf = 0;
	mpz_init(alice.ids.x);
	mpz_init(alice.ids.y);
	alice.ids.inf = 0;
	// --- Anonimização
	mpz_t aux1, aux2;
	mpz_init(aux1);
	mpz_init(aux2);
	mpz_urandomm(alice.k_sess, state, ec.p);
	// gmp_printf("u: %Zd\n", alice.k_sess);
	// ror = Integer(u) | Integer(IDb[0])
	mpz_ior(aux1, alice.k_sess, alice.id.x);
	// gmp_printf("ior: %Zd\n", aux1);
	// st = hashlib.sha256(str(ror)).hexdigest(); rhash = Integer(mod(Integer(st, 16), n))
	hash(aux2, aux1);
	modp(aux2, aux2);
	// gmp_printf("rhash: %Zd\n", aux2);
	mult(&alice.ids, aux2, ec.G);
	// gmp_printf("IDc.x: %Zd\nIDc.y: %Zd\n", alice.ids.x, alice.ids.y);
	// --- Geração de certificado
	mpz_init(alice.v);
	mpz_init(alice.V.x);
	mpz_init(alice.V.y);
	alice.V.inf = 0;
	mpz_urandomm(alice.v, state, ec.p);
	// gmp_printf("v: %Zd\n", alice.v);
	mult(&alice.V, alice.v, ec.G);
	// gmp_printf("V.x: %Zd\nV.y: %Zd\n", alice.V.x, alice.V.y);
	// ror = Integer(IDc[0]) | Integer(V[0])
	mpz_ior(aux1, alice.ids.x, alice.V.x);
	// gmp_printf("ior: %Zd\n", aux1);
	//st = hashlib.sha256(str(ror)).hexdigest(); rhash = Integer(mod(Integer(st, 16), n))
	hash(aux2, aux1);
	modn(aux2, aux2);
	// gmp_printf("rhash: %Zd\n", aux2);
	// sa = mod(mod(ka * rhash, n) + v, n)
	mpz_init(alice.sa);
	mpz_mul(aux1, bob.k, aux2); // aux1 = ka * mod(rhash)
	modn(aux1, aux1); // aux1 = mod(ka * mod(rhash))
	mpz_add(aux2, aux1, alice.v); // aux2 = mod(ka * mod(rhash)) + v
	modp(alice.sa, aux2); // sa = mod(mod(ka * mod(rhash)) + v)
	// gmp_printf("sa: %Zd\n", alice.sa);

	// initialize signature values
	mpz_init(alice.r);
	mpz_init(alice.I.x);
	mpz_init(alice.I.y);
	alice.I.inf = 0;

	mpz_clear(aux1);
	mpz_clear(aux2);
}

inline void sign(struct cert *signedc, int size)
{
	mpz_t aux1, aux2;
	mpz_init(aux1);
	mpz_init(aux2);
	struct coord pt0, pt1;
	mpz_init(pt0.x);
	mpz_init(pt0.y);
	pt0.inf = 0;
	mpz_init(pt1.x);
	mpz_init(pt1.y);
	pt1.inf = 0;
	// - Identidade
	eccAdd(&pt0, alice.ids, alice.id);
	mult(&pt1, alice.k, bob.pu);
	eccAdd(&alice.I, pt0, pt1);
	mpz_clear(pt0.x);
	mpz_clear(pt0.y);
	mpz_clear(pt1.x);
	mpz_clear(pt1.y);
	for (int i = 0; i < size; i++) {
		mpz_urandomm(alice.r, state, ec.p);
		// gmp_printf("m: %Zd\n", signedc[i].m);
		// gmp_printf("r: %Zd\n", alice.r);
		mult(&signedc[i].R, alice.r, ec.G);
		// gmp_printf("R.x: %Zd\nR.y: %Zd\n", signedc[i].R.x, signedc[i].R.y);
		// --- Geração de assinatura
		// s = mod(mod(kb - mod(m * R[0], n), n) * rinv, n)
		mpz_mul(aux1, signedc[i].m, signedc[i].R.x); // aux1 = m * R[0]
		modp(aux1, aux1); // aux1 = mod(m * R[0])
		// gmp_printf("mod(m * R[0]): %Zd\n", aux1);
		mpz_sub(aux2, alice.k, aux1); // aux2 = kb - mod(m * R[0])
		// modp(aux2, aux2); // aux2 = mod(kb - mod(m * R[0]))
		// gmp_printf("mod(kb - mod(m * R[0])): %Zd\n", aux2);
		multInv(aux1, alice.r, ec.n); // rinv = inverse_mod(r, n)
		modn(aux1, aux1);
		// gmp_printf("rinv: %Zd\n", aux1);
		mpz_mul(signedc[i].s, aux2, aux1); // s = mod(kb - mod(m * R[0])) * rinv
		modn(signedc[i].s, signedc[i].s);
		// gmp_printf("s: %Zd\n", signedc[i].s);
	}
	mpz_clear(aux1);
	mpz_clear(aux2);
}

inline int check(struct cert *signedc, int size)
{
	int ret;
	struct coord aux1, aux2, aux3;
	mpz_t s0, s1;
	mpz_init(s0);
	mpz_init(s1);
	mpz_init(aux1.x);
	mpz_init(aux1.y);
	aux1.inf = 0;
	mpz_init(aux2.x);
	mpz_init(aux2.y);
	aux2.inf = 0;
	mpz_init(aux3.x);
	mpz_init(aux3.y);
	aux3.inf = 0;
	// Validação
	mult(&aux1, alice.sa, ec.G); // op1 = Integer(sa) * G
	// gmp_printf("op1.x: %Zd\nop1.y: %Zd\n", aux1.x, aux1.y);
	mpz_ior(s0, alice.ids.x, alice.V.x);
	//st = hashlib.sha256(str(ror)).hexdigest(); rhash = Integer(mod(Integer(st, 16), n))
	hash(s1, s0);
	mult(&aux2, s1, bob.pu);
	eccAdd(&aux3, aux2, alice.V);
	// gmp_printf("op2.x: %Zd\nop2.y: %Zd\n", aux3.x, aux3.y);
	if (mpz_cmp(aux1.x, aux3.x) == 0 && mpz_cmp(aux1.y, aux3.y) == 0) {
		ret = 1;
	} else {
		printf("check: %d: Não validado\n", size);
		ret = 0;
	}
	for (int i = 0; i < size; i++) {
		// Corretude
		// gmp_printf("Kb.x: %Zd\nKb.y: %Zd\n", alice.pu.x, alice.pu.y);
		// op3 = (Integer(s) * R) + (Integer(mod(m * (R[0]), n)) * G)
		mult(&aux1, signedc[i].s, signedc[i].R); // aux1 = (Integer(s) * R)
		// gmp_printf("Integer(s) * R: %Zd, %Zd\n", aux1.x, aux1.y);
		mpz_mul(s0, signedc[i].m, signedc[i].R.x); // s0 = m * (R[0])
		modp(s0, s0); // s0 = mod(m * (R[0])
		// gmp_printf("mod(m * (R[0]): %Zd\n", s0);
		mult(&aux2, s0, ec.G); // aux2 = Integer(mod(m * (R[0]), n)) * G
		// gmp_printf("mod(m * (R[0]), n) * G: %Zd, %Zd\n", aux2.x, aux2.y);
		eccAdd(&aux3, aux1, aux2); // aux3 = Integer(s) * R) + (Integer(mod(m * (R[0]), n)) * G
		// gmp_printf("op3.x: %Zd\nop3.y: %Zd\n", aux3.x, aux3.y);
		if (mpz_cmp(alice.pu.x, aux3.x) == 0 && mpz_cmp(alice.pu.y, aux3.y) == 0) {
			ret = ret*1;
		} else {
			printf("check: %d: Certificado %d incorreto\n", size, i);
			ret = ret*0;
		}
	}
	mpz_clear(s0);
	mpz_clear(s1);
	mpz_clear(aux1.x);
	mpz_clear(aux1.y);
	mpz_clear(aux2.x);
	mpz_clear(aux2.y);
	mpz_clear(aux3.x);
	mpz_clear(aux3.y);
	return ret;
}

void test(long int messageSize) {
	double inicial, final;
	int ret;
	FILE *outp = fopen("output.txt", "a");
	// initialize certificate
	struct cert *signedc = malloc(messageSize * sizeof(struct cert));
	for (int i = 0; i < messageSize; i++) {
		mpz_init(signedc[i].R.x);
		mpz_init(signedc[i].R.y);
		signedc[i].R.inf = 0;
		mpz_init(signedc[i].s);
		mpz_init_set(signedc[i].m, mess[i % nPoints]);
	}

	// Assinatura
	inicial = timestamp();
	sign(signedc, messageSize);
	final = timestamp() - inicial;
	fprintf(outp, "%ld - %f - ", messageSize, final);
	// Verificação
	inicial = timestamp();
	ret = check(signedc, messageSize);
	final = timestamp() - inicial;
	fprintf(outp, "%f\n", final);

	for (int i = 0; i < messageSize; i++) {
		// if (ret == 0) {
		// 	printf("teste %ld: verificação %d inválida\n", messageSize, i);
		// }
		mpz_clear(signedc[i].R.x);
		mpz_clear(signedc[i].R.y);
		mpz_clear(signedc[i].s);
		mpz_clear(signedc[i].m);
	}
	fclose(outp);
}

void main()
{
	// Parameter initialization
	mpz_init_set_str(ec.p, "FFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFF", 16);
	mpz_init_set_str(ec.a, "FFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFC", 16);
	mpz_init_set_str(ec.b, "5AC635D8AA3A93E7B3EBBD55769886BC651D06B0CC53B0F63BCE3C3E27D2604B", 16);
	mpz_init_set_str(ec.G.x, "6B17D1F2E12C4247F8BCE6E563A440F277037D812DEB33A0F4A13945D898C296", 16);
	mpz_init_set_str(ec.G.y, "4FE342E2FE1A7F9B8EE7EB4A7C0F9E162BCE33576B315ECECBB6406837BF51F5", 16);
	ec.G.inf = 0;
	mpz_init_set_str(ec.n, "FFFFFFFF00000000FFFFFFFFFFFFFFFFBCE6FAADA7179E84F3B9CAC2FC632551", 16);

	// Random initialization
	gmp_randinit_mt(state);
	mpz_init_set_str(seed, "C49D360886E704936A6678E1139D26B7819F7E90", 16);
	gmp_randseed(state, seed);

	// mpz_init_set_ui(ec.p, 23);
	// mpz_init_set_ui(ec.a, 1);
	// mpz_init_set_ui(ec.b, 4);
	// mpz_init_set_ui(ec.G.x, 1);
	// mpz_init_set_ui(ec.G.y, 12);
	// ec.G.inf = 0;
	// mpz_init_set_ui(ec.p, 3);

	// gmp_printf("ec.p: %Zd\n", ec.p);
	// gmp_printf("ec.a: %Zd\n", ec.a);
	// gmp_printf("ec.b: %Zd\n", ec.b);
	
	findPoints();
	generateKeys();
	generateCertificates();

	// mpz_t teste;
	// struct coord pt;
	// mpz_init_set_str(teste, "115792089210356248762697446949407573529996955224135760342422259061068512044368", 10);
	// mpz_init(pt.x);
	// mpz_init(pt.y);
	// pt.inf = 0;
	// mult(&pt, teste, ec.G);
	// gmp_printf("%Zd\n%Zd\n", pt.x, pt.y);
	// gmp_printf("G:[%Zd,%Zd]\n", ec.G.x, ec.G.y);
	// gmp_printf("PU[Alice]:[%Zd,%Zd]\n", alice.pu.x, alice.pu.y);
	// gmp_printf("priv[Alice]:%Zd\n", alice.k);
	// gmp_printf("PU[Bob]:[%Zd,%Zd]\n", bob.pu.x, bob.pu.y);
	// gmp_printf("priv[Bob]:%Zd\n", bob.k);

	genMessages();

	// mpz_t res, rop;
	// mpz_init(res);
	// mpz_init(rop);
	// multInv(rop, ec.a);
	// gmp_printf("multinv de %Zd = \n%Zd\n", ec.a, rop);
	// hash(res, alice.k);
	// gmp_printf("hash priv alice: %Zd\n", res);

	test(1);
	test(2);
	test(4);
	test(8);
	test(16);
	test(32);
	test(64);
	test(128);
	test(256);
	test(512);
	test(1024);
	test(2048);
	test(4096);
	test(8192);
	test(16384);
	test(32768);
	test(65536);
	test(131072);
	test(262144);
	test(524288);

	mpz_clear(ec.p);
	mpz_clear(ec.a);
	mpz_clear(ec.b);
	mpz_clear(ec.G.x);
	mpz_clear(ec.G.y);
	mpz_clear(ec.n);
}
