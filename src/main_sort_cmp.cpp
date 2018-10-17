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

//产生随机数函数
void produce_random_number_MT19937(u32* random_array, int size)
{
	int i = 0;
	time_t t; //这里是64位的
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
//希尔排序步长序列,提前计算好
int hibbard_seq[29] = { 536870911, 268435455, 134217727, 67108863, 33554431, 16777215, 8388607, 
				4194303, 2097151, 1048575, 524287, 262143, 131071, 65535, 32767, 16383, 8191, 4095, 2047,
				1023, 511, 255, 127, 63, 31, 15, 7, 3, 1 };
int sedgwick_seq[27] = { 603906049, 268386305, 150958081, 67084289, 37730305, 16764929, 9427969, 4188161,
					2354689, 1045505, 587521, 260609, 146305, 64769, 36289, 16001, 8929, 3905, 2161, 929, 
					505, 209, 109, 41, 19, 5, 1 };
int ciura_seq[25] = { 754868335, 335497038, 149109795, 66271020, 29453787, 13090572, 5818032, 
					2585792, 1149241, 510774, 227011, 100894, 44842, 19930, 8858, 3937,
					/*论文中只到1750，前面的数据按h(k) = INT(2.25 * h(k-1))求得*/1750, 701, 301, 132, 57, 23, 10, 4, 1 };
//一趟希尔插入排序，增量为dk
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
//1959年提出的版本，步长为N/2^k
void shell_sort_shell(u32* A, int length)
{
	int k;
	for (k = length / 2; k > 0; k /= 2)
	{
		shell_insert(A, length, k);
	}
}
//1963年提出的版本，步长为2^k-1
void shell_sort_hibbard(u32* A, int length)
{
	int k = 0;
	while (hibbard_seq[k] >= length) k++;  //找到起始下标
	for (; k < 29; k++)
	{
		shell_insert(A, length, hibbard_seq[k]);
	}
}
//1986年提出的版本，步长为9(4^k-1-2^k/2)+1,4^k+1-6*2^(k+1)/2+1
void shell_sort_sedgwick(u32* A, int length)
{
	int k = 0;
	while (sedgwick_seq[k] >= length) k++;  //找到起始下标
	for (; k < 27; k++)
	{
		shell_insert(A, length, sedgwick_seq[k]);
	}
}
//2001年提出的版本(近似)
void shell_sort_ciura(u32* A, int length)
{
	int k = 0;
	while (ciura_seq[k] >= length) k++;  //找到起始下标
	for (; k < 25; k++)
	{
		shell_insert(A, length, ciura_seq[k]);
	}
}
/****************************************/

/*************quick sort****************/
//partition函数
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
//quick sort主函数
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
//merge函数
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
//merge sort主函数
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
//计数排序
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
		count_array[(A[j] >> offset) & k]--; //注意数组下标比对应位置小1
		B[count_array[(A[j] >> offset) & k]] = A[j];
	}
}
//基数排序
void radix_sort(u32* A, int length, int r)
{
	int loop_cnt = 32 / r;
	u32 max_rbit_num = (1 << r) - 1;
	int i = 0, j = 0, offset = 0;

	if (32 % r != 0) loop_cnt++;  //注意，这里丢了排序结果不对

	for (i = 1; i <= loop_cnt; i++)
	{
		//避免每次将临时数组的结果写回A数组，交替使用A和temp_A存储结果，减少交换次数
		if (i % 2)
			counting_sort(A, temp_array, length, max_rbit_num, offset);
		else
			counting_sort(temp_array, A, length, max_rbit_num, offset);
		offset += r;
	}

	//如果结果还存在temp_A中，将数据写回A中
	if (i % 2 == 0)  //注意这个i最后加一了
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

	//获取待排数据规模
	int n = 0;
	while (1)
	{
		while (1)
		{
			printf("请输入1-2*10^8的待排序规模（-1退出）：");
			scanf("%d", &n);
			getchar();
			if (n == -1) 
				break; //结束
			else if (n < 1 || n > MAX_SIZE)
				printf("数据范围不合法，请重新输入：");
			else
				break;
		}

		if (n == -1)
			break;  //结束

		//耗时变量声明设定
		LARGE_INTEGER n_freq;
		LARGE_INTEGER n_begin_time, n_end_time;
		QueryPerformanceFrequency(&n_freq);

		if (n > 100000)
			printf("待排序数据规模过大，不适合用插入排序\n");
		else
		{
			produce_random_number_MT19937(random_array, n);
			QueryPerformanceCounter(&n_begin_time);
			insertion_sort(random_array, n);
			QueryPerformanceCounter(&n_end_time);
			printf("待排序数据规模为%d, InsertionSort耗时 %f us\n", n, (double)(n_end_time.QuadPart - n_begin_time.QuadPart) * 1000000 / (double)n_freq.QuadPart); //us
		}

		produce_random_number_MT19937(random_array, n);
		QueryPerformanceCounter(&n_begin_time);
		shell_sort_shell(random_array, n);
		QueryPerformanceCounter(&n_end_time);
		printf("待排序数据规模为%d, ShellSort-shell耗时 %f us\n", n, (double)(n_end_time.QuadPart - n_begin_time.QuadPart) * 1000000 / (double)n_freq.QuadPart); //us

		produce_random_number_MT19937(random_array, n);
		QueryPerformanceCounter(&n_begin_time);
		shell_sort_hibbard(random_array, n);
		QueryPerformanceCounter(&n_end_time);
		printf("待排序数据规模为%d, ShellSort-hibbard耗时 %f us\n", n, (double)(n_end_time.QuadPart - n_begin_time.QuadPart) * 1000000 / (double)n_freq.QuadPart); //us

		produce_random_number_MT19937(random_array, n);
		QueryPerformanceCounter(&n_begin_time);
		shell_sort_sedgwick(random_array, n);
		QueryPerformanceCounter(&n_end_time);
		printf("待排序数据规模为%d, ShellSort-sedgwick耗时 %f us\n", n, (double)(n_end_time.QuadPart - n_begin_time.QuadPart) * 1000000 / (double)n_freq.QuadPart); //us

		produce_random_number_MT19937(random_array, n);
		QueryPerformanceCounter(&n_begin_time);
		shell_sort_ciura(random_array, n);
		QueryPerformanceCounter(&n_end_time);
		printf("待排序数据规模为%d, ShellSort-ciura耗时 %f us\n", n, (double)(n_end_time.QuadPart - n_begin_time.QuadPart) * 1000000 / (double)n_freq.QuadPart); //us


		produce_random_number_MT19937(random_array, n);
		QueryPerformanceCounter(&n_begin_time);
		quick_sort(random_array, 0, n - 1);
		QueryPerformanceCounter(&n_end_time);
		printf("待排序数据规模为%d, QuickSort耗时 %f us\n", n, (double)(n_end_time.QuadPart - n_begin_time.QuadPart) * 1000000 / (double)n_freq.QuadPart); //us

		produce_random_number_MT19937(random_array, n);
		QueryPerformanceCounter(&n_begin_time);
		merge_sort(random_array, 0, n - 1);
		QueryPerformanceCounter(&n_end_time);
		printf("待排序数据规模为%d, MergeSort耗时 %f us\n", n, (double)(n_end_time.QuadPart - n_begin_time.QuadPart) * 1000000 / (double)n_freq.QuadPart); //us

		produce_random_number_MT19937(random_array, n);
		QueryPerformanceCounter(&n_begin_time);
		radix_sort(random_array, n, 8);
		QueryPerformanceCounter(&n_end_time);
		printf("待排序数据规模为%d, RadixSort耗时 %f us\n", n, (double)(n_end_time.QuadPart - n_begin_time.QuadPart) * 1000000 / (double)n_freq.QuadPart); //us
	}

#ifdef FILE_TEST
	fclose(fin);
	fclose(fout);
#endif

	return 0;
}