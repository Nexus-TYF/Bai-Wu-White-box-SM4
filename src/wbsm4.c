#include "wbsm4.h"
#include "sbox.h"

M32 L_matrix = {
    .M[0] = 0xA0202080, 
    .M[1] = 0x50101040, 
    .M[2] = 0x28080820, 
    .M[3] = 0x14040410,
    .M[4] = 0xA020208, 
    .M[5] = 0x5010104, 
    .M[6] = 0x2808082, 
    .M[7] = 0x1404041, 
    .M[8] = 0x80A02020, 
    .M[9] = 0x40501010, 
    .M[10] = 0x20280808, 
    .M[11] = 0x10140404, 
    .M[12] = 0x80A0202, 
    .M[13] = 0x4050101, 
    .M[14] = 0x82028080, 
    .M[15] = 0x41014040, 
    .M[16] = 0x2080A020, 
    .M[17] = 0x10405010, 
    .M[18] = 0x8202808, 
    .M[19] = 0x4101404, 
    .M[20] = 0x2080A02, 
    .M[21] = 0x1040501, 
    .M[22] = 0x80820280, 
    .M[23] = 0x40410140, 
    .M[24] = 0x202080A0, 
    .M[25] = 0x10104050, 
    .M[26] = 0x8082028, 
    .M[27] = 0x4041014, 
    .M[28] = 0x202080A, 
    .M[29] = 0x1010405, 
    .M[30] = 0x80808202, 
    .M[31] = 0x40404101
};

void printstate(unsigned char * in)
{
    int i;
    for(i = 0; i < 16; i++) 
    {
        printf("%.2X", in[i]);
    }
    printf("\n");
}

void wbsm4_gen(uint8_t *key)
{
    int i, j, r, x, y;
    uint8_t temp_u8_x, temp_u8_y, temp_u8;
    uint32_t temp_u32;
    uint32_t TD0_u32[6], TD1_u32[6], TD2_u32[3], TR_u32[4], Lc[36], Ec0[32], Ec1[32];
    M32 L[36];
    M32 L_inv[36];
    M8 E[32][2][4];
    M8 E_inv[32][2][4];
    M32 Ei_inv[32][2];
    M32 M[32][2][3];
    M32 C[32];
    M32 LL;

    sm4_context ctx;
    sm4_setkey_enc(&ctx, key);
    InitRandom(((unsigned int)time(NULL)));

    for(r = 0; r < 36; r++) 
    {
        genMatpairM32(&L[r], &L_inv[r]);
        Lc[r] = cus_random();
    }

    for(r = 0; r < 32; r++) 
    {
        for(j = 0; j < 4; j++) 
        {
            genMatpairM8(&E[r][0][j], &E_inv[r][0][j]);
            genMatpairM8(&E[r][1][j], &E_inv[r][1][j]);
        }

        MatrixcomM8to32(E_inv[r][0][0], E_inv[r][0][1], E_inv[r][0][2], E_inv[r][0][3], &Ei_inv[r][0]);
        MatrixcomM8to32(E_inv[r][1][0], E_inv[r][1][1], E_inv[r][1][2], E_inv[r][1][3], &Ei_inv[r][1]);

        MatMulMatM32(Ei_inv[r][0], L_inv[r + 1], &M[r][0][0]);
        MatMulMatM32(Ei_inv[r][0], L_inv[r + 2], &M[r][0][1]);
        MatMulMatM32(Ei_inv[r][0], L_inv[r + 3], &M[r][0][2]);

        MatMulMatM32(Ei_inv[r][1], L_inv[r + 1], &M[r][1][0]);
        MatMulMatM32(Ei_inv[r][1], L_inv[r + 2], &M[r][1][1]);
        MatMulMatM32(Ei_inv[r][1], L_inv[r + 3], &M[r][1][2]);

        MatMulMatM32(L[r + 4], L_inv[r], &C[r]);
    }

    for (r = 0; r < 32; r++)
    {
        MatMulMatM32(L[r + 4], L_matrix, &LL);

        for(i = 0; i < 6; i++)
        {
            TD0_u32[i] = cus_random();
            TD1_u32[i] = cus_random();
        }
        for(i = 0; i < 4; i++)
        {
            TR_u32[i] = cus_random();
        }
         for(i = 0; i < 3; i++)
        {
            TD2_u32[i] = cus_random();
        }
        Ec0[r] = TD0_u32[0] ^ TD0_u32[1] ^ TD0_u32[2] ^ TD0_u32[3] ^ TD0_u32[4] ^ TD0_u32[5];
        Ec1[r] = TD1_u32[0] ^ TD1_u32[1] ^ TD1_u32[2] ^ TD1_u32[3] ^ TD1_u32[4] ^ TD1_u32[5];

        for(x = 0; x < 256; x++) 
        {
            for(j = 0; j < 4; j++)
            {
                temp_u8 = x ^ ((Lc[r] >> (24 - j * 8)) & 0xff);
                temp_u32 = temp_u8 << (24 - j * 8);
                TD[r][0][j][x] = MatMulNumM32(C[r], temp_u32);
            }
            for(j = 0; j < 3; j++)
            {
                TD[r][0][j][x] ^= TD2_u32[j];
            }
            TD[r][0][3][x] ^= Lc[r + 4] ^ TD2_u32[0] ^ TD2_u32[1] ^ TD2_u32[2] ^ TR_u32[0] ^ TR_u32[1] ^ TR_u32[2] ^ TR_u32[3]; 

            for(i = 1; i < 4; i++) 
            {
                temp_u8 = x ^ ((Lc[r + i] >> 24) & 0xff);
                temp_u32 = temp_u8 << 24;
                TD[r][i][0][x] = MatMulNumM32(M[r][0][i - 1], temp_u32);
                temp_u8 = x ^ ((Lc[r + i] >> 16) & 0xff);
                temp_u32 = temp_u8 << 16;
                TD[r][i][1][x] = MatMulNumM32(M[r][0][i - 1], temp_u32);

                temp_u8 = x ^ ((Lc[r + i] >> 8) & 0xff);
                temp_u32 = temp_u8 << 8;
                TD[r][i][2][x] = MatMulNumM32(M[r][1][i - 1], temp_u32);
                temp_u8 = x ^ (Lc[r + i] & 0xff);
                temp_u32 = temp_u8;
                TD[r][i][3][x] = MatMulNumM32(M[r][1][i - 1], temp_u32);
            }

            j = 0;
            for(i = 1; i < 4; i++)
            {
                TD[r][i][0][x] ^= TD0_u32[j++];
                TD[r][i][1][x] ^= TD0_u32[j++];
            }

            j = 0;
            for(i = 1; i < 4; i++)
            {
                TD[r][i][2][x] ^= TD1_u32[j++];
                TD[r][i][3][x] ^= TD1_u32[j++];
            }
        }

        for(x = 0; x < 256; x++)
        {
            for(y = 0; y < 256; y++)
            {
                for(j = 0; j < 4; j++)
                {
                    temp_u8_x = x ^ ((Ec0[r] >> (24 - j * 8)) & 0xff);
                    temp_u8_x = MatMulNumM8(E[r][0][j], temp_u8_x);

                    temp_u8_y = y ^ ((Ec1[r] >> (24 - j * 8)) & 0xff);
                    temp_u8_y = MatMulNumM8(E[r][1][j], temp_u8_y);
                    temp_u8 = SBOX[temp_u8_x ^ temp_u8_y ^ ((ctx.sk[r] >> (24 - j * 8)) & 0xff)];
                    temp_u32 = temp_u8 << (24 - j * 8);
                    TR[r][j][x][y] = MatMulNumM32(LL, temp_u32);
                    TR[r][j][x][y] ^= TR_u32[j];
                }
            }
        }
    }

    //external encoding
    for(i = 0; i < 4; i++) 
    {
        SE[i].Mat = L[i];
        SE[i].Vec.V = Lc[i];

        FE[i].Mat = L_inv[35 - i];
        FE[i].Vec.V = MatMulNumM32(L_inv[35 - i], Lc[35 - i]);
    }
}

void wbsm4_encrypt(unsigned char IN[], unsigned char OUT[])
{
    int r, i, j;
    uint32_t x[36];
    uint32_t s0, s1;
    
    x[0] = GET32(IN);
    x[1] = GET32(IN + 4);
    x[2] = GET32(IN + 8);
    x[3] = GET32(IN + 12);

    for(r = 0; r < 32; r++)
    {
        x[r + 4] = 0;
        s0 = 0;
        s1 = 0;

        for(i = 1; i < 4; i++)
        {
            s0 ^= TD[r][i][0][(x[r + i] >> 24) & 0xff];
            s0 ^= TD[r][i][1][(x[r + i] >> 16) & 0xff];
            s1 ^= TD[r][i][2][(x[r + i] >> 8) & 0xff];
            s1 ^= TD[r][i][3][x[r + i] & 0xff];
        }
        for(j = 0; j < 4; j++)
        {
            x[r + 4] ^= TR[r][j][(s0 >> (24 - j * 8)) & 0xff][(s1 >> (24 - j * 8)) & 0xff];
            x[r + 4] ^= TD[r][0][j][(x[r] >> (24 - j * 8)) & 0xff];
        }
    }

    PUT32(x[35], OUT);
    PUT32(x[34], OUT + 4);
    PUT32(x[33], OUT + 8);
    PUT32(x[32], OUT + 12);
}