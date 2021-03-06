# Sort Algorithm Compare
Sort Algorithm Compare: Insertion Sort, Shell Sort, Quick Sort, Merge Sort and Radix Sort

## 五种排序算法比较

比较Insertion Sort（插入排序），Shell Sort（希尔排序），Quick Sort（快速排序），Merge Sort（归并排序）以及Radix Sort（基数排序）对32位无符号整数的排序效果。

### 1. 数据描述
- 输入数据随机产生，数据范围为 0-2^32-1；
- 输入数据量分别为：10,10^2,10^3,10^4,10^5,10^6,10^7,10^8,2*10^8,10^9。

### 2. 实验环境
- 操作系统：windows10 专业版（64位）
- 处理器：Inter(R) Core(TM) i5-4690K CPU @3.5GHz
- 内存：16G（15.9G可用）
- 编程环境：Visual Studio 2013
- 编程语言：C++
- 编译版本：Release

### 3. 算法分析
插入排序在数据量比较小的时候效率较高，空间开销很小，但平均时间复杂度高，不适用于数据量大的情况；  

希尔排序是直接插入排序的一种改进算法，减少了其复制的次数，速度要快很多，其时间复杂度与步长选取有关【2】； 

快速排序的空间开销也很小，由于其cache命中率很高，使其与现代计算机的体系结构有很好的契合性，效率很高；  

归并排序采用分治思想，将待排序数据分为左右两部分，递归排序，再将排序结果合并到一个集合中，其需要一个大小为n的临时存储空间用以保存合并序列；  

基数排序与前面四种排序算法不同，并不是基于比较的排序算法，它对待排序数的每一位进行计数排序，计数排序的运行时间为O(n)，在一定条件下，基数排序的效率很高，可达O(n)。

五种算法的特性如下表所示。

 ![image](https://github.com/chenwhql/SortCompare/blob/master/images/sortfeature.png)
 

### 4. 算法实现
#### 4.1  随机数生成
本实验要求生成随机数的数据范围为0~2^32-1，而C++的rand()函数只能产生0~RAND_MAX(如在Visual Studio 2013上面最大为0x7FFF，即32767)范围内的随机数，生成 范围内的随机数常用的有三种方法【1】：

（1） Mersenne twister生成算法：号称是目前最好的随机数生成算法，它是由Takuji Nishimura 和 Makoto Matsumoto于1997年开发的一种随机数生成方法，它基于有限二进制字段上的矩阵线性再生，可以快速产生高质量的伪随机数，该算法的循环周期为2^19937-1。

（2） Windows API CryptGenRandom方法：这个函数是windows用来提供给生成随机密码数字功能所用，一般应用在随机给vector赋值或者生成噪声点时使用。这个函数所提供的随机性要好于rand()函数。

（3） rand()函数拼接法：使用rand()函数随机生成的三个数来拼接成一个32位数字的方法，这种方法实现起来也很简单，但是由于rand()的随机性并不是很好，所以不推荐这种方法。一个32位随机数由(rand()<<17)|(rand()<<2)|(rand())方法拼合而成。

本实验采用方法（1）Mersenne twister生成算法的变种Mersenne Twister MT19937算法，保证测试数据集有较好的随机性。

#### 4.2  Insertion Sort实现
按照《算法导论》翻译版第10页的伪码实现，不再赘述。

#### 4.3  Shell Sort实现
希尔排序实现时选取不同步长序列，其算法效率也不同，如表所示【2】。

 ![image](https://github.com/chenwhql/SortCompare/blob/master/images/shellsort.png)

 
本实验中选取四种比较有代表性的步长序列进行算法性能比较。


#### 4.4  Quick Sort实现
按照《算法导论》翻译版第95页的伪码实现，不再赘述。

#### 4.5  Merge Sort实现
基本按照《算法导论》翻译版第17页的伪码实现，区别有两点：
（1） 在Merge时以左右数组的长度为结束标志，而不是将最后一个值设为无穷；
（2） 若左右任意一侧元素全部被处理，则结束比较，将另一侧数据全部加入尾部。

#### 4.6  Radix Sort实现
按照《算法导论》翻译版第109页与111页的伪码实现，主要思路是将32位二进制数看作几个r位数拼接而成，然后分别对每一段r位二进制数进行计数排序。  
根据《算法导论》的描述，直接使用理论值欠妥，首先需要对这个结论进行验证，确认基数排序中最佳的r取值。  
这里在实现的时候，为了减少过程中数据复制的次数，使用目标数组A与临时数组B交替存放过程排序结果。  

### 5. 实验过程
### 5.1  数据采集
因为实现时采用的time(&t)作为种子，每次生成的随机数都是不用的。所以在取r值和排序的耗时数据之前，对于每一个数量级，都先生成对应数量级的随机数写入到文件中，目的是在每次排序后可以从文件中读取并恢复数组的初始状态，以保证初始条件的一致性。

每次针对一个数量级进行测试，生成数据后存入excle表格。
### 5.2  最终可执行程序
最终实现的可执行程序目的是便于展示，与采集数据时所用的程序有以下区别：
（1） 最终可执行程序在调用每个排序算法前，随机数集都是新生成的；
（2） 最终可执行程序没有文件读写操作。
最终实现效果如图所示。

 ![image](https://github.com/chenwhql/SortCompare/blob/master/images/result.png)
 

### 6. 实验结果及分析
### 6.1  Radix Sort算法r值分析
取 ，分别测试其在不同输入规模下计数排序的耗时，计算平均数据整理如表所示。表中同一输入数量级下使用的数据样本相同，每个r值对应的耗时值均是三次实验数据的平均值（见附件数据分析.xlsx），单位为us，具有较高的准确性与可对比性。

 ![image](https://github.com/chenwhql/SortCompare/blob/master/images/radixsort_rval.png)
 
将表中数据绘制成折线图，如图所示。
 
 ![image](https://github.com/chenwhql/SortCompare/blob/master/images/radixsort_rval_line.png)

从图表中可以看出，在n规模较小时，r的最佳取值出现在4-6区间，但当输入规模达到 时，r的最佳取值便收敛到8，不再浮动，总体而言，r最佳实验值为8，此时基数排序的耗时最短，效率最高。下表是理论值与实验值的对比。

 ![image](https://github.com/chenwhql/SortCompare/blob/master/images/rval_comp.png)
 
由表可知，当n规模很大时，r的理论取值便不再适用于实际情况。其主要原因可能是受限于计算机高速cache的大小，当r值较大时，使得计数排序中k的取值增大，辅助数组的空间增大，一次对r位数的排列便需要占用较大的cache空间，当达到一定限度时，cache的空间不能满足一次计数排序所需要的辅助空间大小，排序的效率便开始下降。
综上，这里选择r=8作为最终算法使用值。

### 6.2  排序算法比较
分别测试共计8种排序算法在不同输入规模下的排序耗时，计算平均数据整理如下表所示。表中同一输入数量级下使用的数据样本相同，每个r值对应的耗时值均是三次实验数据的平均值（见附件数据分析.xlsx），单位为us。这里插入排序只计算到 ， 运行时间过长，单方面拉高纵坐标轴，不利于折线图的观察，所以输入规模大于等于 之后不再考虑插入排序。

 ![image](https://github.com/chenwhql/SortCompare/blob/master/images/timecost.png)
 
将上表数据绘制成折线图，如图所示。

 ![image](https://github.com/chenwhql/SortCompare/blob/master/images/sort_compare.png)

分析数据，可以总结出以下观点：
（1） n较小时（ 规模以下），插入排序性能最好，希尔排序（sedgwick步长序列与ciura步长序列）也有很好的表现。
（2） n到达 规模以上时，插入排序的效果直线下降， 规模以上便不再适用，同时基数排序与快速排序的优势开始变得明显。


（3） 希尔排序选择4种不同的步长，当数据输入规模足够大时，效率差距明显，使用sedgwick序列与ciura序列的希尔排序算法耗时大约只有使用普通步长序列时的一半，实验证明 序列效率并没有比 序列高多少，甚至在数据规模较大时，表现更差。
（4） 快速排序和归并排序虽然理论上有着相同的时间复杂度，相同的平均效率，但是实际来看，快速排序的效果要比归并排序好很多，关键在于归并排序空间复杂度较高，有较多的赋值操作，其cache命中率不如快速排序。
（5） 基数排序非比较的排序方式显示出很好的效果，即使在数据达到 规模时，依旧耗时很短，其所耗时间比快速排序小一个数量级。
（6） n数量级足够大时，按时间效率从高到低对算法进行排序，序列如下：基数排序>快速排序>归并排序>希尔排序>插入排序。

### 7. 实验结论
1. 基数排序r值得选取要考虑到硬件的限制，理论值不一定是最佳选择，像本实验中的最佳值为8，并不等于lgn；
2. 插入排序小规模数据排序效率很高，不适合于大规模数据的排序；
3. 希尔排序采用不同的步长序列，会有差异较为明显的表现，在数据规模适中时，采用合适序列的希尔排序效率可与快速排序媲美，但是数据规模较大时，依旧不如快速排序；
4. 快速排序很适合于大规模数据的排序，相比同样时间复杂度的归并排序有着更好的表现，但是如果在空间条件允许的情况下，基数排序的表现优于快速排序，因此条件允许的情况下基数排序也是很好的选择。

### 参考资料
[1] 使用C++如何产生32位随机数，http://www.cnblogs.com/thu539/archive/2011/11/14/2247717.html
[2] Shellsort, https://en.wikipedia.org/wiki/Shellsort

