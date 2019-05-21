/*
 *实现论文 “Erasure Coding in Windows Azure Storage” 里面介绍的Local Reconstruction Codes算法
 *Author：黄曙权（shughuangqqq@gmail.com)
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "LRC.h"


#define dmalloc(type, num) (type *) malloc(sizeof(type)*(num))
#define mult(x,y) mult_table[x][y]
#define div(x,y) div_table[x][y]


 /*
  *这个结构存储具体的每次操作
  *如果ope == 0 那么就把 data[s][si]的值复制到 data[d][di]的位置
  *如果ope == 1 那么 data[d][di] ^= data[s][si]
  */
typedef struct {
	unsigned char ope, s, si, d, di;
} operation;


/*
 *具体的编码和解码操作
 *need_para需要的参数，从0~lrc_n 依次表示数据块，组编码块，全局编码块， 值为1表示需要，0表示不需要
 *tot_ope 为ope数据的长度
 *ope为一个数组，存储每一步需要的操作。
 */
typedef struct {
	int *need_para;
	int tot_ope;
	operation *ope;
} task;


/*
 *encode_task 存储编码对应的task, 其以所需的源数据块的顺序进行排序。
 *decode_task_cache 存储每一个不同情况的解码task。
 *每一种需要decode都可以根据state变成一个数值（详见lrc_state_to_long(int*))
 *如果decode[type] == NULL那么表示type对应的数据损坏情况不能解码。
 *这样就decode_task_cache的长度为2^lrc_n次方，而其中非NULL的个数大约为C(lrc_n,lrc_l+lrc_r);
 *如果lrc_n比较大的时候，可以改变type的表示方法，从而减少长度，不过只要lrc_n在20（大约需要8M内存）以内，一般是可以接受的大小
 *encode_offset 存储每一个数据块在编码操作对应到encode_task里面ope位置。 也就是encode_offset[i] 到 encode_offset[i+1] -1 对应的ope会用到数据块i
 */
static task *encode_task;
static task **decode_task_cache;
static int encode_offset[lrc_k + 1];

/*
 *temp_var 对应程序初始化需要的一些变量，lrc_init()最后的free_temp_var()就会释放这些内存。
 *mult_table 存储GF(2^8)的乘法表
 *div_table 存储GF(2^8)的除法表
 *encode_matrix 存储初始的编码矩阵
 *bitmap[k]存储一个为二维矩阵，表示k这个数在GF(2^8)中对应的bitmap
 *tot_bit[k]表示k这个数据对应的bitmap中有多少个1（方便统计个数和申请空间）
 *block_task[k]表示一个k数据块编码时的task。 所以encode_task 就是所有block_task的合并
 *now_ope对应每个数据中的bit位，如果为0，表示这个时候应该复制，否则就是做xor运算
 */
typedef struct {
	gf **mult_table;
	gf **div_table;
	int *encode_matrix;
	unsigned char ***bitmap;
	unsigned char *tot_bit;
	task *block_task;
	int **now_ope;
} temp_var;




/*
 *生成GF(2^8)的乘法和除法表
 *本源多项式为x^8 + x^4 + x^3 + x^2 + 1  （这个多项式直接来自 http://blog.csdn.net/mengboy/article/details/1514445 这篇博客）
 *a是为了计算方便初始化a^i对应的值
 */
static int make_mult_table(temp_var *var) {
	gf a[] = { 1,2,4,8,16,32,64,128,29,58,116,232,205,135,19 };
	int i, j, k;
	gf ** mult_table = dmalloc(gf*, 256);
	if (mult_table == NULL) {
		return LRC_ERROR_OUT_OF_MEMORY;
	}
	gf ** div_table = dmalloc(gf*, 256);
	if (div_table == NULL) {
		free(mult_table);
		return LRC_ERROR_OUT_OF_MEMORY;
	}
	memset(mult_table, 0, sizeof(gf*) * 256);
	memset(div_table, 0, sizeof(gf*) * 256);
	for (i = 0; i < 256; i++) {
		mult_table[i] = dmalloc(gf, 256);
		div_table[i] = dmalloc(gf, 256);
		if (mult_table[i] == NULL || div_table[i] == NULL) {
			for (; i >= 0; i--) {
				if (mult_table[i] != NULL) free(mult_table[i]);
				if (div_table[i] != NULL) free(div_table[i]);
			}
			free(mult_table);
			free(div_table);
			return LRC_ERROR_OUT_OF_MEMORY;
		}

	}
	unsigned int z;
	for (i = 0; i < 256; i++)
		for (j = 0; j < 256; j++) {
			z = 0;
			//多项式乘法
			for (k = 0; k < 8; k++)
				if (j & (1 << k))  z ^= (i << k);
			mult_table[i][j] = 0;
			//多项式乘法结果取余(低位保留，大于2**8的高位化解取余）
			for (k = 0; k < 15; k++)
				mult_table[i][j] ^= ((z & (1 << k)) == 0) ? 0 : a[k];
			//检测交换律结果
			if (i > j && mult_table[i][j] != mult_table[j][i]) {
				fprintf(stderr, "error: %d %d %d %d\n", i, j, mult_table[i][j], mult_table[j][i]);
				exit(-1);
			}

			//从乘法表中得到除法表
			z = mult_table[i][j];
			if (i != 0 && j != 0) {
				div_table[z][i] = j;
				div_table[z][j] = i;
			}
		}
	var->mult_table = mult_table;
	var->div_table = div_table;
	return 0;
}


/*
 *把now_ope这个数组设置为0；
 */
static void reset_now_ope(temp_var *var) {
	int i, j;
	for (i = 0; i < lrc_n; i++)
		for (j = 0; j < 8; j++)
			var->now_ope[i][j] = 0;
}


/*
 *把data[d] = data[d] xor (k * data[s]) 的运算变成bitmap下的运算序列，并存储在now里面
 *在bitmap的情况下 data[d][i] = sum(foreach(j) (k[i][j]*data[s][j])), GF(2^8)里面，sum 操作等同于xor操作
 *因为k[i][j] 只有0和1两种情况，且当为0是0*data[s][j] == 0  ==> data[d][i]不变，所以0的运算可以忽略
 *而k[i][j]为1 的时候 k[i][j]*data[s][j] == data[s][j] 这样就不需要乘法，直接计算sum, 也就是xor。
 *又因为一开始data[d][i]需要赋值为0，那么就相当第一次碰到k[i][j] == 1时，data[d][i] = data[s][j]; 这是ope为0，而此后所有运算应是xor运算
 */
static operation* schedule_to_operation(temp_var *var, operation* now, int s, int d, int k) {
	if (k == 0) return now;
	int x, y;
	for (x = 0; x < 8; x++) {
		for (y = 0; y < 8; y++) {
			if (var->bitmap[k][x][y] == 1) {
				now->s = s;
				now->si = y;
				now->d = d;
				now->di = x;
				now->ope = var->now_ope[d][x];
				var->now_ope[d][x] = 1;
				now++;
			}
		}
	}
	return now;
}


/*
 *首先构造根据MS 论文里的方法构造出编码矩阵
 *然后把矩阵转换为一个操作串
 *同时构造了block_task
 */
static int encoding_matrix(temp_var *var) {
	int i, j, rs, k;
	int *encode_matrix = dmalloc(int, lrc_k * lrc_n);
	if (encode_matrix == NULL) return LRC_ERROR_OUT_OF_MEMORY;
	var->encode_matrix = encode_matrix;
	gf **mult_table = var->mult_table;
	gf **div_table = var->div_table;
	for (i = 0; i < lrc_k; i++) {
		for (j = 0; j < lrc_k; j++)
			encode_matrix[i*lrc_k + j] = 0;
		encode_matrix[i*lrc_k + i] = 1;
	}
	for (i = lrc_k; i < lrc_k + lrc_l; i++) {
		for (j = 0; j < lrc_k; j++)
			encode_matrix[i*lrc_k + j] = ((i - lrc_k)*lrc_group_size <= j && j < (i - lrc_k + 1)*lrc_group_size) ? 1 : 0;
	}
	int *temp = dmalloc(int, lrc_k);
	if (temp == NULL) return LRC_ERROR_OUT_OF_MEMORY;
	int *l;
	l = (encode_matrix + (lrc_k + lrc_l)*lrc_k);
	for (i = 0; i < lrc_l; i++) {
		for (j = 1; j <= lrc_group_size; j++) {
			temp[i*lrc_group_size + j - 1] = j << (i * 4);
			*(l++) = *(temp + i * lrc_group_size + j - 1);
		}
	}
	for (i = 1; i < lrc_r; i++) {
		for (j = 0; j < lrc_k; j++) {
			*l = mult(temp[j], *(l - lrc_k));
			l++;
		}
	}
	free(temp);

	// matrix to task;
	encode_task = dmalloc(task, 1);
	if (encode_task == NULL) return LRC_ERROR_OUT_OF_MEMORY;
	memset(encode_task, 0, sizeof(task));
	task *block_task = dmalloc(task, lrc_n);
	if (block_task == NULL) return LRC_ERROR_OUT_OF_MEMORY;
	memset(block_task, 0, lrc_n * sizeof(task));
	var->block_task = block_task;
	for (i = 0; i < lrc_n; i++) block_task[i].ope = NULL;
	unsigned char *tot_bit = var->tot_bit;
	encode_task->need_para = dmalloc(int, lrc_n);
	if (encode_task->need_para == NULL) return LRC_ERROR_OUT_OF_MEMORY;
	for (i = 0; i < lrc_n; i++) encode_task->need_para[i] = i < lrc_k ? 1 : 0;
	encode_task->tot_ope = 0;
	reset_now_ope(var);
	operation *now;

	for (i = lrc_k; i < lrc_n; i++) {
		rs = i * lrc_k;
		block_task[i].tot_ope = 0;
		for (j = 0; j < lrc_k; j++) {
			if (encode_matrix[rs + j] != 0)
				block_task[i].tot_ope += tot_bit[encode_matrix[rs + j]];
		}
		block_task[i].ope = dmalloc(operation, block_task[i].tot_ope);
		if (block_task[i].ope == NULL) return LRC_ERROR_OUT_OF_MEMORY;
		now = block_task[i].ope;
		for (j = 0; j < lrc_k; j++) {
			if (encode_matrix[rs + j] != 0)
				now = schedule_to_operation(var, now, j, i, encode_matrix[rs + j]);
		}
		encode_task->tot_ope += block_task[i].tot_ope;
	}
	encode_task->ope = dmalloc(operation, encode_task->tot_ope);
	if (encode_task->ope == NULL) return LRC_ERROR_OUT_OF_MEMORY;
	now = encode_task->ope;
	for (i = lrc_k; i < lrc_n; i++) {
		memcpy(now, block_task[i].ope, sizeof(operation)*block_task[i].tot_ope);
		now += block_task[i].tot_ope;
	}
	//利用ope的源地址进行排序，并产生对应的encode_offset, 以便于对于一个连续的src block进行多次操作，增加CPU cache locality
	//由于tot_ope并不大。就直接进行选择排序
	for (i = 0; i < lrc_k; i++) encode_offset[i] = -1;
	for (i = 0; i < encode_task->tot_ope; i++) {
		k = i;
		for (j = i + 1; j < encode_task->tot_ope; j++)
			if (encode_task->ope[j].s < encode_task->ope[k].s || (encode_task->ope[j].s == encode_task->ope[k].s && encode_task->ope[j].si < encode_task->ope[k].si))
				k = j;
		if (i != k) {
			operation temp;
			temp = encode_task->ope[i];
			encode_task->ope[i] = encode_task->ope[k];
			encode_task->ope[k] = temp;
		}
		if (encode_offset[encode_task->ope[i].s] == -1) encode_offset[encode_task->ope[i].s] = i;
	}
	encode_offset[lrc_k] = encode_task->tot_ope;

	return 0;
}


/*
 *矩阵求逆。结果存放到inv里面
 *算法：
 *    对要求逆的矩阵A 产生矩阵 AE （E为单位矩阵）然后根据行间的运算，转化为 EB 矩阵，那么矩阵B即为A的逆
 */
static int invert_matrix(temp_var *var, int *mat, int *inv, int rows) {
	int i, j, k, x, rs2, rs, temp, inverse;
	gf **mult_table = var->mult_table;
	gf **div_table = var->div_table;
	k = 0;
	//初始化为E。
	for (i = 0; i < rows; i++) {
		for (j = 0; j < rows; j++) {
			inv[k] = (i == j) ? 1 : 0;
			k++;
		}
	}

	//对每一行进行操作
	for (i = 0; i < rows; i++) {
		rs = rows * i;

		//如果mat[i][i] 为0，那就找出mat[j][i]不为0的j，并加到i行。 如果找不到，那说明该矩阵不能求逆。退出
		if (mat[rs + i] == 0) {
			for (j = i + 1; j < rows && mat[rows*j + i] == 0; j++);
			if (j == rows) return -1;
			rs2 = j * rows;
			for (k = 0; k < rows; k++) {
				mat[rs + k] ^= mat[rs2 + k];
				inv[rs + k] ^= inv[rs2 + k];
			}
		}

		k = rs + i;
		//把mat[i][i]变为1；
		temp = mat[k];
		if (temp != 1) {
			inverse = div(1, temp);
			for (j = 0; j < rows; j++) {
				mat[rs + j] = mult(mat[rs + j], inverse);
				inv[rs + j] = mult(inv[rs + j], inverse);
			}
		}

		//消去mat[j][i] (j!=i)的值
		k = -rows + i;
		for (j = 0; j < rows; j++) {
			k += rows;
			if (j == i) continue;
			if (mat[k] != 0) {
				if (mat[k] == 1) {
					rs2 = rows * j;
					for (x = 0; x < rows; x++) {
						mat[rs2 + x] ^= mat[rs + x];
						inv[rs2 + x] ^= inv[rs + x];
					}
				}
				else {
					temp = mat[k];
					rs2 = rows * j;
					for (x = 0; x < rows; x++) {
						mat[rs2 + x] ^= mult(temp, mat[rs + x]);
						inv[rs2 + x] ^= mult(temp, inv[rs + x]);
					}
				}
			}

		}

	}
	return 0;
}


/*
 *根据提供的参数构造矩阵，生成解码矩阵
 */
static int decoding_matrix(temp_var *var, int *inv, int rows, int *selected_rows, int *selected_cols) {
	int *encode_matrix = var->encode_matrix;
	int *mat = dmalloc(int, rows*rows);
	if (mat == NULL) return LRC_ERROR_OUT_OF_MEMORY;
	int i, j, rs, rs2;
	for (i = 0; i < rows; i++) {
		rs = i * rows;
		rs2 = selected_rows[i] * lrc_k;
		for (j = 0; j < rows; j++) mat[rs + j] = encode_matrix[rs2 + selected_cols[j]];
	}
	int error = invert_matrix(var, mat, inv, rows);
	if (error) {
		free(mat);
		return error;
	}
	free(mat);
	return 0;
}


/*
 *根据提供的选择行和列的信息，构建新解码矩阵，并生成对应的操作串
 *need_decode需要解码的标记信息，其实就是损坏块的state
 */
static int construct_operation(temp_var *var, int rows, int *selected_rows, int *selected_cols, int *need_decode, task* t) {
	int i, j, k, rc;
	int *inv = dmalloc(int, rows*rows);
	if (inv == NULL) return LRC_ERROR_OUT_OF_MEMORY;

	int error = decoding_matrix(var, inv, rows, selected_rows, selected_cols);
	if (error) {
		free(inv);
		return LRC_ERROR_OUT_OF_MEMORY;
	}

	//根据编码矩阵产生操作串
	reset_now_ope(var);
	unsigned char *tot_bit = var->tot_bit;
	t->tot_ope = 0;
	for (i = 0; i < rows; i++) {
		if (need_decode[selected_cols[i]] == 0) continue;
		rc = i * rows;
		for (j = 0; j < rows; j++) if (inv[rc + j] != 0) t->tot_ope += (tot_bit[inv[rc + j]]);
	}
	t->ope = dmalloc(operation, t->tot_ope);
	if (t->ope == NULL) {
		free(inv);
		return LRC_ERROR_OUT_OF_MEMORY;
	}
	operation *now = t->ope;
	for (i = 0; i < rows; i++) {
		if (need_decode[selected_cols[i]] == 0) continue;
		rc = i * rows;
		for (j = 0; j < rows; j++)
			if (inv[rc + j] != 0) {
				now = schedule_to_operation(var, now, selected_rows[j], selected_cols[i], inv[rc + j]);
			}
	}
	free(inv);
	return 0;
}

/*
 *根据state，算出需要的参数。并判断能否解码。
 *如果能解码就构建对应的task
 *need_para 需要的数据块标记。
 *selected_group 那些小组被选择的标记
 */
static int state_to_task(temp_var *var, int *state, int seat, int tot) {
	int i, j, k, temp;
	int error;
	int *need_para = dmalloc(int, lrc_n);
	if (need_para == NULL) return LRC_ERROR_OUT_OF_MEMORY;
	int *selected_group = dmalloc(int, lrc_l);
	if (selected_group == NULL) { free(need_para); return LRC_ERROR_OUT_OF_MEMORY; }
	unsigned char *tot_bit = var->tot_bit;
	task *block_task = var->block_task;
	for (i = 0; i < lrc_n; i++) need_para[i] = 0;

	// 以下找出需要的处理的数据组，和统计需要全局编码块的个数    
	k = 0;
	temp = 0;  // 统计需要多少个全局编码块
	for (i = 0; i < lrc_l; i++) {
		int b = 0;  // b相当于一个bool类型，如果这个数据组中有损坏数据，那么就b为true，否者为false
		for (j = 0; j < lrc_group_size; j++) b |= state[k + j];
		b |= (state[lrc_k + i]);
		selected_group[i] = b;  // 如果有数据快损坏，那么表示这个数据组是需要被处理的。
		if (b) {
			for (j = 0; j < lrc_group_size; j++)
				if (state[k + j] == 0) {      //如果这个数据没有损坏，那么这个数据就需要被提供来解码
					need_para[k + j] = 1;
				}
				else {
					if (need_para[lrc_k + i] == 0 && state[lrc_k + i] == 0) {   //如果小组编码块没有损坏，并且没有被标记为需要数据块，那么就可以用小组编码块来解码顺坏数据块
						need_para[lrc_k + i] = 1;
					}
					else {   //否则需要一个全局数据块才可以解码
						temp++;
					}
				}
		}
		k += lrc_group_size;
	}
	for (i = lrc_k + lrc_l; i < lrc_n; i++) {   //如果有一个全局编码块被顺坏，那当需要一个全局数据块
		if (state[i] == 1) temp++;
	}

	// 如果需要的全局编码块个数超过全局编码块个数，不能被解码。直接退出。
	if (temp > lrc_r) {
		free(need_para);
		free(selected_group);
		return 0;
	}

	// 如果需要全局数据块参与，那就意味着所有小组的数据块都需要提供来解码。
	if (temp > 0) {
		for (i = 0; i < lrc_l; i++)
			if (selected_group[i] == 0) {
				selected_group[i] = 1;
				k = i * lrc_group_size;
				for (j = 0; j < lrc_group_size; j++) need_para[k + j] = 1;
			}
	}
	for (i = lrc_k + lrc_l; i < lrc_n; i++) {
		if (state[i] == 1) temp--;
	}
	i = lrc_k + lrc_l;
	while (temp > 0) {
		if (state[i] == 0) {
			need_para[i] = 1;
			temp--;
		}
		i++;
	}

	//申请需要的空间
	int type = lrc_state_to_long(state);
	task *t = dmalloc(task, 1);
	if (t == NULL) {
		free(need_para);
		free(selected_group);
		return LRC_ERROR_OUT_OF_MEMORY;
	}
	t->need_para = need_para;
	int rows = 0;
	for (i = 0; i < lrc_l; i++) rows += ((selected_group[i] == 0) ? 0 : 1)*lrc_group_size;
	int *selected_rows = dmalloc(int, rows);
	if (selected_rows == NULL) {
		free(need_para);
		free(selected_group);
		free(t);
		return LRC_ERROR_OUT_OF_MEMORY;
	}
	int *selected_cols = dmalloc(int, rows);
	if (selected_cols == NULL) {
		free(need_para);
		free(selected_group);
		free(t);
		free(selected_rows);
		return LRC_ERROR_OUT_OF_MEMORY;
	}

	// 计算需要的行和列
	k = 0;
	for (i = 0; i < lrc_n; i++)
		if (need_para[i] == 1) {
			selected_rows[k++] = i;
		}
	k = 0;
	for (i = 0; i < lrc_l; i++)
		if (selected_group[i] != 0) {
			int l = i * lrc_group_size;
			for (j = 0; j < lrc_group_size; j++) {
				selected_cols[k++] = l + j;
			}
		}

	error = construct_operation(var, rows, selected_rows, selected_cols, state, t);
	if (error) {
		free(need_para);
		free(selected_group);
		free(t);
		free(selected_rows);
		free(selected_cols);
		return LRC_ERROR_OUT_OF_MEMORY;
	}

	temp = 0; //多少个编码块损坏
	for (i = lrc_k; i < lrc_n; i++) {
		if (state[i] == 1) {
			temp += block_task[i].tot_ope;
		}
	}

	//如果有编码块损坏，那么只要直接读取之前预处理的block_task,合并到刚才计算的task里面，就可以完成解码。
	if (temp > 0) {
		operation *tope = dmalloc(operation, t->tot_ope + temp);
		if (tope == NULL) {
			free(need_para);
			free(selected_group);
			free(t->ope);
			free(t);
			free(selected_rows);
			free(selected_cols);
			return LRC_ERROR_OUT_OF_MEMORY;
		}
		operation *now = tope;
		memcpy(tope, t->ope, t->tot_ope * sizeof(operation));
		now += t->tot_ope;
		free(t->ope);
		t->ope = tope;
		t->tot_ope += temp;
		for (i = lrc_k; i < lrc_n; i++) {
			if (state[i] == 1) {
				memcpy(now, block_task[i].ope, sizeof(operation)*block_task[i].tot_ope);
				now += block_task[i].tot_ope;
			}
		}
	}
	decode_task_cache[type] = t;
	free(selected_group);
	free(selected_rows);
	free(selected_cols);
	return 0;
}


/*
 *搜索数据块损坏的情况。 使用深搜
 *state为一个长度lrc_n的数组，1表示对应的数据块已经坏了
 *seat表示目前搜索到第几块
 *tot表示目前一共有多少块是损坏的。 tot不能超过lrc_l+lrc_r
 */
static int search_state(temp_var *var, int *state, int seat, int tot) {
	if (seat >= lrc_n) return 0;
	if (tot == lrc_r + lrc_l) return 0;
	int error;
	error = search_state(var, state, seat + 1, tot);
	if (error) return error;
	state[seat] = 1;
	tot++;
	error = search_state(var, state, seat + 1, tot);
	if (error) return error;
	error = state_to_task(var, state, seat, tot);
	if (error) return error;
	state[seat] = 0;
	return 0;
}



/*
 *构建decode_task_cache
 */
static int construct_cache(temp_var *var) {
	int i = 0;
	decode_task_cache = dmalloc(task*, (1 << lrc_n));
	if (decode_task_cache == NULL) return LRC_ERROR_OUT_OF_MEMORY;
	for (i = 0; i < (1 << lrc_n); i++) {
		decode_task_cache[i] = NULL;
	}
	int *state = dmalloc(int, lrc_n);
	if (state == NULL) return LRC_ERROR_OUT_OF_MEMORY;
	for (i = 0; i < lrc_n; i++) state[i] = 0;
	int error = search_state(var, state, 0, 0);
	if (error) {
		free(state);
		return error;
	}
	free(state);
	return 0;
}


/*
 *初始化bitmap,tot_bit,和申请now_ope的空间
 *k转换成对应的bitmap，相关算法参考Jerasure。
 *具体算法： bitmap[k][i] 对应k*2^i的二进制位
 */
static int init_bitmap(temp_var *var) {
	gf **mult_table = var->mult_table;
	gf **div_table = var->div_table;

	unsigned char ***bitmap = dmalloc(unsigned char **, 256);
	if (bitmap == NULL) return LRC_ERROR_OUT_OF_MEMORY;
	memset(bitmap, 0, 256 * sizeof(unsigned char **));
	var->bitmap = bitmap;

	unsigned char *tot_bit = dmalloc(unsigned char, 256);
	if (tot_bit == NULL) return LRC_ERROR_OUT_OF_MEMORY;
	memset(tot_bit, 0, 256 * sizeof(unsigned char));
	var->tot_bit = tot_bit;

	int **now_ope = dmalloc(int *, lrc_n);
	if (now_ope == NULL) return LRC_ERROR_OUT_OF_MEMORY;
	var->now_ope = now_ope;
	memset(now_ope, 0, lrc_n * sizeof(int*));

	int i, j, k, x, y;
	for (i = 0; i < 256; i++) {
		tot_bit[i] = 0;
		gf k = i;
		bitmap[i] = dmalloc(unsigned char *, 8);
		if (bitmap[i] == NULL) return LRC_ERROR_OUT_OF_MEMORY;
		memset(bitmap[i], 0, sizeof(unsigned char *) * 8);
		for (x = 0; x < 8; x++) {
			bitmap[i][x] = dmalloc(unsigned char, 8);
			if (bitmap[i][x] == NULL) return LRC_ERROR_OUT_OF_MEMORY;
			for (y = 0; y < 8; y++) {
				bitmap[i][x][y] = ((k &(1 << y)) ? 1 : 0);
				tot_bit[i] += bitmap[i][x][y];
			}
			k = mult(k, 2);
		}
	}

	for (i = 0; i < lrc_n; i++) {
		now_ope[i] = dmalloc(int, 8);
		if (now_ope[i] == NULL) return LRC_ERROR_OUT_OF_MEMORY;
		for (j = 0; j < 8; j++) now_ope[i][j] = 0;
	}
	return 0;
}


/*
 *释放var的内存。 需要遍历一遍var里面所有指针，保证都被释放
 */
static void free_temp_var(temp_var *var) {
	if (var != NULL) {
		int i, x;
		if (var->mult_table != NULL) {
			for (i = 0; i < 256; i++)
				if (var->mult_table[i] != NULL) free(var->mult_table[i]);
			free(var->mult_table);
		}
		if (var->div_table != NULL) {
			for (i = 0; i < 256; i++)
				if (var->div_table[i] != NULL) free(var->div_table[i]);
			free(var->div_table);
		}
		if (var->encode_matrix != NULL) free(var->encode_matrix);
		if (var->tot_bit != NULL) free(var->tot_bit);
		if (var->bitmap != NULL) {
			for (i = 0; i < 256; i++)
				if (var->bitmap[i] != NULL) {
					for (x = 0; x < 8; x++) if (var->bitmap[i][x] != NULL) free(var->bitmap[i][x]);
					free(var->bitmap[i]);
				}
			free(var->bitmap);
		}
		if (var->block_task != NULL) {
			for (i = 0; i < lrc_n; i++)
				if (var->block_task[i].ope != NULL) free(var->block_task[i].ope);
			free(var->block_task);
		}
		if (var->now_ope != NULL) {
			for (i = 0; i < lrc_n; i++) {
				if (var->now_ope[i] != NULL) free(var->now_ope[i]);
			}
			free(var->now_ope);
		}
		free(var);
	}
}

/*
 *根据给定的operation序列，对数据进行编解码。
 *本程序速度的核心。
 *以packet为单位进行编解码。
 *对于一个packet，均分为8块，设为b0,b2,...,b7。分别对应这个数据的bit位的数据。
 *(b0i,b1i,...,b7i)对应packet的第i个数据 （b0i表示b0第i个bit位）。
 * 根据矩阵的并行性，就可以进行加速运算。每次取出sizeof（long）个数据一起运算。（也就是cpu字长）
 * 由于数据方法的切割以packet为单位，而不是每个byte进行，因此不同的packet_size会导致不同的计算结果。
 */
static void bit_code(gf **data, task *t, long size) {
	int i, j;
	long *s, *d, *end;
	operation *ope;
	for (i = 0; i < size; i += lrc_packet_size) {
		ope = t->ope;
		for (j = 0; j < t->tot_ope; j++) {
			s = (long *)(data[ope->s] + i + (ope->si << lrc_bitpacket_skip));
			d = (long *)(data[ope->d] + i + (ope->di << lrc_bitpacket_skip));
			if (ope->ope == 1) {
				end = (long *)((char *)s + lrc_bitpacket_size);
				while (s < end) {
					(*d) ^= (*s);
					s++;
					d++;
				}
			}
			else {
				memcpy(d, s, lrc_bitpacket_size);
			}
			ope++;
		}
	}
}





int lrc_init() {
	lrc_free();              // 避免之前内存没有释放，而造成内存泄露。

	temp_var *var = dmalloc(temp_var, 1);
	if (var == NULL) return LRC_ERROR_OUT_OF_MEMORY;
	memset(var, 0, sizeof(temp_var));
	int error = make_mult_table(var);
	if (error) {
		free_temp_var(var);
		lrc_free();
		return error;
	}
	error = init_bitmap(var);
	if (error) {
		free_temp_var(var);
		lrc_free();
		return error;
	}
	error = encoding_matrix(var);
	if (error) {
		free_temp_var(var);
		lrc_free();
		return error;
	}
	error = construct_cache(var);
	if (error) {
		free_temp_var(var);
		lrc_free();
		return error;
	}
	free_temp_var(var);
	return 0;
}

void lrc_free() {
	int i;
	if (encode_task != NULL) {
		if (encode_task->need_para != NULL) free(encode_task->need_para);
		if (encode_task->ope != NULL) free(encode_task->ope);
		free(encode_task);
	}
	if (decode_task_cache != NULL) {
		for (i = 0; i < (1 << lrc_n); i++)
			if (decode_task_cache[i] != NULL) {
				if (decode_task_cache[i]->need_para != NULL) free(decode_task_cache[i]->need_para);
				if (decode_task_cache[i]->ope != NULL) free(decode_task_cache[i]->ope);
				free(decode_task_cache[i]);
			}
		free(decode_task_cache);
	}
}


void lrc_encoding(gf **data, long size) {
	bit_code(data, encode_task, size);
}

void lrc_encoding_part(gf **newData, gf **oldData, long size) {
	int i, j, k;
	long *s, *d, *end, *olds;
	operation *ope, *endope;
	for (i = 0; i < size; i += lrc_packet_size) {
		//Copy the old coding blocks
		for (j = lrc_k; j < lrc_n; j++)
			memcpy(newData[j] + i, oldData[j] + i, lrc_packet_size);

		for (k = 0; k < lrc_k; k++)
			if (newData[k] != NULL) {
				ope = encode_task->ope + encode_offset[k];
				endope = encode_task->ope + encode_offset[k + 1];
				for (; ope < endope; ope++) {
					//no matter ope is 1(xor) or 0(copy), we should offset the impact of old data block,
					//and add the impact of new data block. Since the add and decrease in GF() is xor, so
					//the final result is:  code_block = code_block ^ old_data_block ^ new_data_block
					s = (long *)(newData[ope->s] + i + (ope->si << lrc_bitpacket_skip));
					olds = (long *)(oldData[ope->s] + i + (ope->si << lrc_bitpacket_skip));
					d = (long *)(newData[ope->d] + i + (ope->di << lrc_bitpacket_skip));
					end = (long *)((char *)s + lrc_bitpacket_size);
					while (s < end) {
						(*d) ^= (*s) ^ (*olds);
						s++;
						d++;
						olds++;
					}

				}
			}
	}

}


long lrc_state_to_long(int *state) {
	int i;
	long type = 0;
	for (i = 0; i < lrc_n; i++) type |= (long)state[i] << i;
	return type;
}

int *lrc_need_para(long type) {
	if (decode_task_cache[type] == NULL) return NULL;
	else return decode_task_cache[type]->need_para;
}

void lrc_decoding(long type, gf **data, long size) {
	bit_code(data, decode_task_cache[type], size);
}











