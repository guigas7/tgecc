#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sha256.c"

/*
Output should be:
ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad
248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1
cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0
*/

void print_hash(unsigned char hash[]) {
    int idx;
    for (idx = 0; idx < 32; idx++)
        printf("%02x", hash[idx]);
    printf("\n");
}

int main() {
    unsigned char
        text2[] = {
            "FFFFFFFF00000001000000000000000000000000FFFFFFFFFFFFFFFFFFFFFFFC"
        },
        hash[32];
    int size2 = 256;
    int res;

    SHA256_CTX ctx;

    // Hash two
    sha256_init(&ctx);
    sha256_update(&ctx, text2, strlen(text2));
    sha256_final(&ctx, hash, &res, size2, 23);
    printf("hash: \n");
    print_hash(hash);

    // getchar();
    return 0;
}