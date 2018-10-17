#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <windows.h>
#include "mt19937ar.h"
using namespace std;

//#define LOCAL
//#define FILE_TEST
#define MAX_SIZE 200000000
#define CNT_SIZE 65536

typedef unsigned int u32;

u32 random_array[MAX_SIZE];
u32 temp_array[MAX_SIZE];
u32 count_array[CNT_SIZE];

//�������������
void produce_random_number_MT19937(u32* random_array, int size)
{
	int i = 0;
	time_t t; //������64λ��
	init_genrand((u32)time(&t));
	for (i = 0; i < size; i++)
	{
		random_array[i] = genrand_int32();
	}
}

/**********insertion sort**************/
void insertion_sort(u32* A, int length)
{
	int i, j, key;
	for (j = 1; j < length; j++)
	{
		key = A[j];
		i = j - 1;
		while (i >= 0 && A[i] > key)
		{
			A[i + 1] = A[i];
			i = i - 1;
		}
		A[i + 1] = key;
	}
}
/**************************************/

/************shell sort****************/
//ϣ�����򲽳�����,��ǰ�����
int hibbard_seq[29] = { 536870911, 268435455, 134217727, 67108863, 33554431, 16777215, 8388607, 
				4194303, 2097151, 1048575, 524287, 262143, 131071, 65535, 32767, 16383, 8191, 4095, 2047,
				1023, 511, 255, 127, 63, 31, 15, 7, 3, 1 };
int sedgwick_seq[27] = { 603906049, 268386305, 150958081, 67084289, 37730305, 16764929, 9427969, 4188161,
					2354689, 1045505, 587521, 260609, 146305, 64769, 36289, 16001, 8929, 3905, 2161, 929, 
					505, 209, 109, 41, 19, 5, 1 };
int ciura_seq[25] = { 754868335, 335497038, 149109795, 66271020, 29453787, 13090572, 5818032, 
					2585792, 1149241, 510774, 227011, 100894, 44842, 19930, 8858, 3937,
					/*������ֻ��1750��ǰ������ݰ�h(k) = INT(2.25 * h(k-1))���*/1750, 701, 301, 132, 57, 23, 10, 4, 1 };
//һ��ϣ��������������Ϊdk
void shell_insert(u32* A, int length, int dk)
{
	int i, j, key;
	for (j = dk; j < length; j++)
	{
		key = A[j];
		i = j - dk;
		while (i >= 0 && A[i] > key)
		{
			A[i + dk] = A[i];
			i = i - dk;
		}
		A[i + dk] = key;
	}
}
//1959������İ汾������ΪN/2^k
void shell_sort_shell(u32* A, int length)
{
	int k;
	for (k = length / 2; k > 0; k /= 2)
	{
		shell_insert(A, length, k);
	}
}
//1963������İ汾������Ϊ2^k-1
void shell_sort_hibbard(u32* A, int length)
{
	int k = 0;
	while (hibbard_seq[k] >= length) k++;  //�ҵ���ʼ�±�
	for (; k < 29; k++)
	{
		shell_insert(A, length, hibbard_seq[k]);
	}
}
//1986������İ汾������Ϊ9(4^k-1-2^k/2)+1,4^k+1-6*2^(k+1)/2+1
void shell_sort_sedgwick(u32* A, int length)
{
	int k = 0;
	while (sedgwick_seq[k] >= length) k++;  //�ҵ���ʼ�±�
	for (; k < 27; k++)
	{
		shell_insert(A, length, sedgwick_seq[k]);
	}
}
//2001������İ汾(����)
void shell_sort_ciura(u32* A, int length)
{
	int k = 0;
	while (ciura_seq[k] >= length) k++;  //�ҵ���ʼ�±�
	for (; k < 25; k++)
	{
		shell_insert(A, length, ciura_seq[k]);
	}
}
/****************************************/

/*************quick sort****************/
//partition����
u32 partition(u32* A, int p, int r)
{
	int i, j;
	u32 x = A[r];
	i = p - 1;
	for (j = p; j < r; j++)
	{
		if (A[j] <= x)
		{
			i = i + 1;
			swap(A[i], A[j]);
		}
	}
	swap(A[i + 1], A[r]);
	return i + 1;
}
//quick sort������
void quick_sort(u32* A, int p, int r)
{
	int q = 0;
	if (p < r)
	{
		q = partition(A, p, r);
		quick_sort(A, p, q-1);
		quick_sort(A, q+1, r);
	}
}
/****************************************/

/************merge sort******************/
//merge����
void merge(u32* A, int p, int q, int r)
{
	int i, j, k;
	for (i = p; i <= q; i++)
	{
		temp_array[i] = A[i];
	}
	for (j = q + 1; j <= r; j++)
	{
		temp_array[j] = A[j];
	}
	i = p;
	j = q + 1;
	k = p;
	while (i <= q && j <= r)
	{
		if (temp_array[i] <= temp_array[j])
		{
			A[k++] = temp_array[i++];
		}
		else
		{
			A[k++] = temp_array[j++];
		}
	}

	while (i <= q) A[k++] = temp_array[i++];
	while (j <= r) A[k++] = temp_array[j++];
}
//merge sort������
void merge_sort(u32* A, int p, int r)
{
	int q = 0;
	if (p < r)
	{
		q = (p + r) / 2;
		merge_sort(A, p, q);
		merge_sort(A, q + 1, r);
		merge(A, p, q, r);
	}
}
/****************************************/

/***************radix sort***************/
//��������
void counting_sort(u32* A, u32* B, int length, u32 k, u32 offset)
{
	int i, j;
	for (i = 0; i <= k; i++)
	{
		count_array[i] = 0;
	}
	for (j = 0; j < length; j++)
	{
		count_array[(A[j] >> offset) & k]++;
	}
	for (i = 1; i <= k; i++)
	{
		count_array[i] = count_array[i] + count_array[i - 1];
	}
	for (j = length - 1; j >= 0; j--)
	{
		count_array[(A[j] >> offset) & k]--; //ע�������±�ȶ�Ӧλ��С1
		B[count_array[(A[j] >> offset) & k]] = A[j];
	}
}
//��������
void radix_sort(u32* A, int length, int r)
{
	int loop_cnt = 32 / r;
	u32 max_rbit_num = (1 << r) - 1;
	int i = 0, j = 0, offset = 0;

	if (32 % r != 0) loop_cnt++;  //ע�⣬���ﶪ������������

	for (i = 1; i <= loop_cnt; i++)
	{
		//����ÿ�ν���ʱ����Ľ��д��A���飬����ʹ��A��temp_A�洢��������ٽ�������
		if (i % 2)
			counting_sort(A, temp_array, length, max_rbit_num, offset);
		else
			counting_sort(temp_array, A, length, max_rbit_num, offset);
		offset += r;
	}

	//������������temp_A�У�������д��A��
	if (i % 2 == 0)  //ע�����i����һ��
	{
		for (j = 0; j < length; j++)
		{
			A[j] = temp_array[j];
		}
	}
}
/****************************************/


int main()
{
#ifdef LOCAL
	freopen("10.txt","r",stdin);
	freopen("out.txt","w",stdout);
#endif

#ifdef FILE_TEST
	FILE *fin, *fout;
	fin = fopen("10_000_000.txt","rb");
	fout = fopen("out.txt","wb");
#endif

	//��ȡ�������ݹ�ģ
	int n = 0;
	while (1)
	{
		while (1)
		{
			printf("������1-2*10^8�Ĵ������ģ��-1�˳�����");
			scanf("%d", &n);
			getchar();
			if (n == -1) 
				break; //����
			else if (n < 1 || n > MAX_SIZE)
				printf("���ݷ�Χ���Ϸ������������룺");
			else
				break;
		}

		if (n == -1)
			break;  //����

		//��ʱ���������趨
		LARGE_INTEGER n_freq;
		LARGE_INTEGER n_begin_time, n_end_time;
		QueryPerformanceFrequency(&n_freq);

		if (n > 100000)
			printf("���������ݹ�ģ���󣬲��ʺ��ò�������\n");
		else
		{
			produce_random_number_MT19937(random_array, n);
			QueryPerformanceCounter(&n_begin_time);
			insertion_sort(random_array, n);
			QueryPerformanceCounter(&n_end_time);
			printf("���������ݹ�ģΪ%d, InsertionSort��ʱ %f us\n", n, (double)(n_end_time.QuadPart - n_begin_time.QuadPart) * 1000000 / (double)n_freq.QuadPart); //us
		}

		produce_random_number_MT19937(random_array, n);
		QueryPerformanceCounter(&n_begin_time);
		shell_sort_shell(random_array, n);
		QueryPerformanceCounter(&n_end_time);
		printf("���������ݹ�ģΪ%d, ShellSort-shell��ʱ %f us\n", n, (double)(n_end_time.QuadPart - n_begin_time.QuadPart) * 1000000 / (double)n_freq.QuadPart); //us

		produce_random_number_MT19937(random_array, n);
		QueryPerformanceCounter(&n_begin_time);
		shell_sort_hibbard(random_array, n);
		QueryPerformanceCounter(&n_end_time);
		printf("���������ݹ�ģΪ%d, ShellSort-hibbard��ʱ %f us\n", n, (double)(n_end_time.QuadPart - n_begin_time.QuadPart) * 1000000 / (double)n_freq.QuadPart); //us

		produce_random_number_MT19937(random_array, n);
		QueryPerformanceCounter(&n_begin_time);
		shell_sort_sedgwick(random_array, n);
		QueryPerformanceCounter(&n_end_time);
		printf("���������ݹ�ģΪ%d, ShellSort-sedgwick��ʱ %f us\n", n, (double)(n_end_time.QuadPart - n_begin_time.QuadPart) * 1000000 / (double)n_freq.QuadPart); //us

		produce_random_number_MT19937(random_array, n);
		QueryPerformanceCounter(&n_begin_time);
		shell_sort_ciura(random_array, n);
		QueryPerformanceCounter(&n_end_time);
		printf("���������ݹ�ģΪ%d, ShellSort-ciura��ʱ %f us\n", n, (double)(n_end_time.QuadPart - n_begin_time.QuadPart) * 1000000 / (double)n_freq.QuadPart); //us


		produce_random_number_MT19937(random_array, n);
		QueryPerformanceCounter(&n_begin_time);
		quick_sort(random_array, 0, n - 1);
		QueryPerformanceCounter(&n_end_time);
		printf("���������ݹ�ģΪ%d, QuickSort��ʱ %f us\n", n, (double)(n_end_time.QuadPart - n_begin_time.QuadPart) * 1000000 / (double)n_freq.QuadPart); //us

		produce_random_number_MT19937(random_array, n);
		QueryPerformanceCounter(&n_begin_time);
		merge_sort(random_array, 0, n - 1);
		QueryPerformanceCounter(&n_end_time);
		printf("���������ݹ�ģΪ%d, MergeSort��ʱ %f us\n", n, (double)(n_end_time.QuadPart - n_begin_time.QuadPart) * 1000000 / (double)n_freq.QuadPart); //us

		produce_random_number_MT19937(random_array, n);
		QueryPerformanceCounter(&n_begin_time);
		radix_sort(random_array, n, 8);
		QueryPerformanceCounter(&n_end_time);
		printf("���������ݹ�ģΪ%d, RadixSort��ʱ %f us\n", n, (double)(n_end_time.QuadPart - n_begin_time.QuadPart) * 1000000 / (double)n_freq.QuadPart); //us
	}

#ifdef FILE_TEST
	fclose(fin);
	fclose(fout);
#endif

	return 0;
}