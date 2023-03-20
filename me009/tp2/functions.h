void remplirDirichlet(float An[], float As[], float Ap[], float Ae[], float Aw[], float b[],
                      int nn, int nx, int ny,
                      float gs, float ge, float gw, float gn)
{
    for (int i = 0; i < nn; i++)
    {
        Ap[i] = 4;
        if (i < (ny - 1) * nx)
            An[i] = -1;
        if (i >= nx)
            As[i] = -1;
        if (((i + 1) % nx) != 0)
            Ae[i] = -1;
        if (((i) % nx) != 0)
            Aw[i] = -1;

        if (i < nx)
            b[i] += gs;
        if (i >= (ny - 1) * nx)
            b[i] += gn;
        if ((i % nx) == 0)
            b[i] += gw;
        if (((i + 1) % nx) == 0)
            b[i] += ge;
    }
}
void initQ(float Q[], int nn)
{
    for (int i = 0; i < nn; i++)
    {
        Q[i] = 1;
    }
}
void bMinusAq(float An[], float As[], float Ap[], float Ae[], float Aw[], float b[],
              int nn, int nx, int ny, float Q[], float r[])
{
    for (int i = 0; i < nn; i++)
    {
        r[i] = b[i] - Ap[i] * Q[i];
        if (i > 0)
            r[i] -= Aw[i] * Q[i - 1];
        if (i > nx - 1)
            r[i] -= As[i] * Q[i - nx];
        if (i < (ny - 1) * nx)
            r[i] -= An[i] * Q[i + nx];
        if (i < (ny)*nx - 1)
            r[i] -= Ae[i] * Q[i + 1];
    }
    // printf("r:\n");
    // for (int i = 0; i < nn; i++)
    //     printf("%.4f, ", r[i]);
    // printf("\n");
}
void descentedq(float An[], float As[], float Ap[], float Ae[], float Aw[], float b[],
                int nn, int nx, int ny, float r[], float dq[], float omega)
{
    for (int i = 0; i < nn; i++)
    {
        dq[i] = omega * r[i];
        if (i > 0)
            dq[i] -= omega * Aw[i] * dq[i - 1];
        if (i >= nx)
            dq[i] -= omega * As[i] * dq[i - nx];
        dq[i] = dq[i] / Ap[i];
    }
//     printf("dq:\n");
//     for (int i = 0; i < nn; i++)
//         printf("%.4f, ", dq[i]);
//     printf("\n");
}
//   (An, As, Ap, Ae, Aw, b, nn, nx, ny, r, dQ, omega);

void amendQ(float dQ[], float Q[], int nn)
{
    for (int i = 0; i < nn; i++)
    {
        Q[i] = Q[i] + dQ[i];
    }
}