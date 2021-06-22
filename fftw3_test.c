#include<stdio.h>
#include<math.h>
#include <time.h>
#include<stdlib.h>
#include<fftw3.h>
#include<complex.h>
#define N 4096*105
//typedef struct{
//        double real;
//        double imag;
//}complex;
float complex x[N], *W;
#define size 4096
double l[100000],h[100000];
double PI=(4.0*atan(1));




float check_time(char flag)
{
    static clock_t start;
    clock_t temp;
    float ret_time = 0;
    temp = clock();
    if(flag != 0){
        ret_time = (float)(temp - start)/CLOCKS_PER_SEC;
    }
    start = temp;
    return ret_time;
}


int load_hex_file(char *filename,char *data_buf,int data_len)
{
        FILE * fp;
        fp = fopen(filename,"r");
        if(fp == NULL)
        {
                printf("open file error\n");
                exit(-1);
        }

        fgets(data_buf,data_len,fp);
        fclose(fp);
}

int save_file(char *filename,char *data_buf,int data_len)
{
        FILE * fp;
        fp = fopen(filename,"r");
        if(fp == NULL)
        {
                printf("open file error\n");
                exit(-1);
        }

        fgets(data_buf,data_len,fp);
        fclose(fp);
}


void  change()
{
        float complex temp;
        unsigned short i=0,j=0,k=0;
        double t;
        for(i=0;i<size;i++)
        {
                k=i;
                j=0;
                t=(log(size)/log(2));
                while( (t--)>0 )
                {
                        j=j<<1;
                        j|=(k & 1);
                        k=k>>1;
                }
                if(j>i)
                {
                        temp=x[i];
                        x[i]=x[j];
                        x[j]=temp;
                }
        }
}
void transform()
{
        int i;
        W=(float complex *)malloc(sizeof(float complex) * size);
        for(i=0;i<size;i++)
        {
            W[i] = cos(2*PI/size*i) + -1*sin(2*PI/size*i)*I;
        }
}

void fft()
{
    int i=0,j=0,k=0,m=0;
    float complex q,y,z;
    change();
    for(i=0;i<log(size)/log(2) ;i++)
    {
        m=1<<i;
        for(j=0;j<size;j+=2*m)
        {
            for(k=0;k<m;k++)
            {
                q = x[k+j+m] * W[size*k/2/m];
                y = x[j+k] + q;
                z = x[j+k] + q;
                x[j+k]=y;
                x[j+k+m]=z;
            }
        }
    }
}

void fftw3_test(char * data_buffer,int data_len)
{
    #define ROUND_TIMES 100
    fftw_complex *in, *out;
    fftw_plan p;
    
    double fft_amp[ROUND_TIMES*size];
    
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    for(int i = 0 ;i < data_len;i ++)
    {
        in[i][0] = data_buffer[i*2];
        in[i][1] = data_buffer[i*2 + 1];
//        if(i%4096 == 0)
//            printf("data in is %f + %fi\n",in[i][0],in[i][1]);
    }
    
    check_time(0);
    for(int i = 0;i < ROUND_TIMES;i ++)
    {
        int base_idx = i * size;
        
        p = fftw_plan_dft_1d(4096, &in[base_idx], &out[base_idx], FFTW_FORWARD, FFTW_ESTIMATE);
        
        fftw_execute(p); /* repeat as needed */

        for(int j = 0;j < size;j ++)
        {
            double temp = out[j+base_idx][0]*out[j+base_idx][0] + out[j+base_idx][1]*out[j+base_idx][1];
            fft_amp[j + base_idx] = (sqrt(temp));
            //printf("fft data out %f + %fi\n",);
        }
    } 
    double data_max = 0;
    for(int i = 0;i < data_len;i ++)
    {
        if(data_max < fft_amp[i])
            data_max = fft_amp[i];
    }
  
    for(int i = 0;i < data_len;i ++)
    {
        fft_amp[i] = 20*log10(fft_amp[i]/data_max);
    }  
    printf("fftw3 run time:%lfs\n",check_time(1));
    
    unsigned char fft_int8[ROUND_TIMES*size*3]; 
    unsigned int temp_int;
    double ufft_val;
    FILE * fp;
    fp = fopen("fftw3_result.bin","w");
    for(int i = 0;i < data_len;i ++)
    {
        ufft_val = (fft_amp[i]/data_max * 0xFFFFF);
        temp_int = (unsigned int)ufft_val ;
//        fft_int8[i*3 + 0] = temp_int&0xFF;
//        fft_int8[i*3 + 1] = (temp_int>>8)&0xFF;
//        fft_int8[i*3 + 2] = (temp_int>>16)&0xFF;
        fprintf(fp,"%02f ",fft_amp[i]);
        if(i%size == size-1)
        {
           fprintf(fp,"\n"); 
        }
    }    
    
    fclose(fp);
    
    
    fftw_destroy_plan(p);
    fftw_free(in); fftw_free(out);
}





int main()
{
    int i;
    double fs;
    double t[100000];
    clock_t begin, end;
    double cost1, cost2, persent;
    #define DATA_LEN 4096*2
    #define ROUND_TIMES 100
    char data_buf[DATA_LEN*ROUND_TIMES + 10];
    
    fs = 2.6*1000000;
    load_hex_file("FMcapture1.dat",data_buf,DATA_LEN*100);
    for(int i  = 0;i < DATA_LEN*100;i ++)
    {
        data_buf[i] = data_buf[i] - 127;
    }
    
    fftw3_test(data_buf,ROUND_TIMES*size);
    
    /*****************************/
    for(i=0;i<size;i++)
    {
            t[i]=i/fs;  //时间序列
            h[i]=i*(fs/(double)size);//频率序列
    }
    transform();
    
    check_time(0);
    for(int i = 0;i < ROUND_TIMES;i ++)
    {
        int base_idx = i * size;
        //printf("base_idx[%d] is: %d\n", i,base_idx);
        for(unsigned int j = 0;j < size;j ++)
        {
            unsigned int data_idx = base_idx + j*2;
            x[i + data_idx] = (data_buf[base_idx + data_idx    ] ) + (data_buf[base_idx + data_idx + 1])*I;
        }
        fft();
    }
    //printf("myfft run time:%lfs\n",check_time(1));

    /*fft shift*/
    //check_time(0); 
    for(int i = 0;i < ROUND_TIMES;i ++)
    {
        int base_idx = i * size;
        //printf("base idx is : %d\n",base_idx);
        for(int j = 0;j < size;j ++)
        {
            int idx = base_idx + j;
            float temp;
            double ffreal,ffimg;
            ffreal = creal(x[idx]);
            ffimg = cimag(x[idx]);
            //temp = sqrt(x[idx].imag*x[idx].imag+x[idx].real*x[idx].real);
            temp = sqrt(ffreal*ffreal + ffimg*ffimg);
            if(j < size/2){
                l[base_idx + size/2 + j] = temp;
            }else{
                l[base_idx - size/2 + j] = temp;
            }
        }
    }
    
    double data_max = 0;
    for(int i = 0;i < ROUND_TIMES*size; i++)
    {
        if(data_max < l[i])
            data_max = l[i];
        
    }
    
    for(int i = 0;i < ROUND_TIMES*size; i++)
    {
        l[i] = 20*log10(l[i]/data_max);
        
    }

    printf("run time:%lfs\n",check_time(1));
    FILE * fp;
    fp = fopen("myfft_result.bin","w");
    for(int i = 0;i < ROUND_TIMES*size; i++)
    {
        double temp;
        temp = l[i];
        fprintf(fp,"%02f ",temp);
        if(i%size == size-1)
        {
           fprintf(fp,"\n"); 
        }        
    }
    
    fclose(fp);
    
    


    
    return 0;
}



