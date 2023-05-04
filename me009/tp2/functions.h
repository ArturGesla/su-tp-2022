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
void remplirNeumann(float An[], float As[], float Ap[], float Ae[], float Aw[], float b[],
                    int nn, int nx, int ny,
                    float gs, float ge, float gw, float gn, float h)
{
    gw = 0;
    ge = 1;
    gs = 0;
    gn = 0;
    for (int i = nx; i < nn - nx; i++)
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
        {
            b[i] += ge * ((i + 1) / nx - 1) * h;
            printf("%d \t %lf \n", i, b[i]);
        }
    }

    // tob and bottom stripe to fill

    for (int i = 0; i < nx; i++)
    {
        Ap[i] = 2;
        if (i < (ny - 1) * nx)
            An[i] = -1;
        if (i >= nx)
            As[i] = -1;
        if (((i + 1) % nx) != 0)
            Ae[i] = -1 / 2.0;
        if (((i) % nx) != 0)
            Aw[i] = -1 / 2.0;

        if (i < nx)
            b[i] += gs;
        if (i >= (ny - 1) * nx)
            b[i] += gn;
        if ((i % nx) == 0)
            b[i] += gw / 2.0;
        if (((i + 1) % nx) == 0)
        {
            b[i] += ge * ((i + 1) / nx - 1) / 2.0 * h;
            printf("%d \t %lf \n", i, b[i]);
        }
    }
    for (int i = nn - nx; i < nn; i++)
    {
        Ap[i] = 2;
        if (i < (ny - 1) * nx)
            An[i] = -1;
        if (i >= nx)
            As[i] = -1;
        if (((i + 1) % nx) != 0)
            Ae[i] = -1 / 2.0;
        if (((i) % nx) != 0)
            Aw[i] = -1 / 2.0;

        if (i < nx)
            b[i] += gs;
        if (i >= (ny - 1) * nx)
            b[i] += gn;
        if ((i % nx) == 0)
            b[i] += gw / 2.0;
        if (((i + 1) % nx) == 0)
        {
            b[i] += ge * ((i + 1) / nx - 1) / 2.0 * h;
            printf("%d \t %lf \n", i, b[i]);
        }
    }
}
float diffWithAnalytical(float Q[], int nn, int nx, int ny, float h, int modes)
{
    float normDiff=0;
    for (int ix = 0; ix < nx; ix++)
    {
        for (int iy = 0; iy < ny; iy++)
        {
            float x = (ix + 1) * h;
            float y = iy * h;
            float l=1.25;
            float L=1.5;
            float analytical = 1/2.0*l/L*x;
            float pi=3.14159265359;
            for (int i = 1; i < modes; i++)
            // {
            //     float bn=2.0/l/sinh(i*pi)*l*l/i/i/pi/pi*(pow(-1,i)-1);
            //     printf("current bn: %4.2e\n",bn );
            //     analytical+=bn*cos(y*i*pi/l)*sinh(i*pi*x/l);
            // }
            {
                float bn=2.0/l*l*l/i/i/pi/pi*(pow(-1,i)-1);
                printf("current bn: %4.2e\n",bn );
                analytical+=bn*cos(y*i*pi/l)*exp(i*pi*(x/l-1));
            }
            //trzebaby tu pomyslec
            Q[ix+iy*nx]=analytical;
                        // Q[ix+iy*nx]=Q[ix+iy*nx]-analytical;

            printf("Diff at ix %d iy %d is %4.2e\n",ix,iy,Q[ix+iy*nx]-analytical);
                     normDiff+=   (Q[ix+iy*nx]-analytical)*(Q[ix+iy*nx]-analytical)/nx/ny;

            
        }
    }
    return sqrt(normDiff);
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
float fluxV4Dirichlet(float Q[], int nx, int ny, float h)
{
    float QV = 0;
    int iy = (ny - 1) / 2;
    printf("iy: %d\n", iy);
    for (int ix = 0; ix < nx + 1; ix++)
    {
        if (ix == 0)
        {
            float dTdyym = 0;
            float dTdyy = (-Q[ix + (iy - 1) * nx] + Q[ix + (iy + 1) * nx]) / 2.0 / h;
            // printf("ix: %d\tdTdyym: %lf\tdTdy: %lf\n", ix, dTdyym, dTdyy);

            QV += (dTdyym + dTdyy) / 2.0 * h;
        }
        else if (ix == nx)
        {
            float dTdyym = (-Q[ix - 1 + (iy - 1) * nx] + Q[ix - 1 + (iy + 1) * nx]) / 2.0 / h;
            float dTdyy = 0;
            // printf("ix: %d\tdTdyym: %lf\tdTdy: %lf\n", ix, dTdyym, dTdyy);

            QV += (dTdyym + dTdyy) / 2.0 * h;
        }
        else
        {
            float dTdyym = (-Q[ix - 1 + (iy - 1) * nx] + Q[ix - 1 + (iy + 1) * nx]) / 2.0 / h;
            float dTdyy = (-Q[ix + (iy - 1) * nx] + Q[ix + (iy + 1) * nx]) / 2.0 / h;
            // printf("ix: %d\tdTdyym: %lf\tdTdy: %lf\n", ix, dTdyym, dTdyy);

            QV += (dTdyym + dTdyy) / 2.0 * h;
        }
    }

    return -QV;
}