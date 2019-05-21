/*
 *ʵ������ ��Erasure Coding in Windows Azure Storage�� ������ܵ�Local Reconstruction Codes�㷨
 *Author������Ȩ��shughuangqqq@gmail.com)
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "LRC.h"


#define dmalloc(type, num) (type *) malloc(sizeof(type)*(num))
#define mult(x,y) mult_table[x][y]
#define div(x,y) div_table[x][y]


 /*
  *����ṹ�洢�����ÿ�β���
  *���ope == 0 ��ô�Ͱ� data[s][si]��ֵ���Ƶ� data[d][di]��λ��
  *���ope == 1 ��ô data[d][di] ^= data[s][si]
  */
typedef struct {
	unsigned char ope, s, si, d, di;
} operation;


/*
 *����ı���ͽ������
 *need_para��Ҫ�Ĳ�������0~lrc_n ���α�ʾ���ݿ飬�����飬ȫ�ֱ���飬 ֵΪ1��ʾ��Ҫ��0��ʾ����Ҫ
 *tot_ope Ϊope���ݵĳ���
 *opeΪһ�����飬�洢ÿһ����Ҫ�Ĳ�����
 */
typedef struct {
	int *need_para;
	int tot_ope;
	operation *ope;
} task;


/*
 *encode_task �洢�����Ӧ��task, ���������Դ���ݿ��˳���������
 *decode_task_cache �洢ÿһ����ͬ����Ľ���task��
 *ÿһ����Ҫdecode�����Ը���state���һ����ֵ�����lrc_state_to_long(int*))
 *���decode[type] == NULL��ô��ʾtype��Ӧ��������������ܽ��롣
 *������decode_task_cache�ĳ���Ϊ2^lrc_n�η��������з�NULL�ĸ�����ԼΪC(lrc_n,lrc_l+lrc_r);
 *���lrc_n�Ƚϴ��ʱ�򣬿��Ըı�type�ı�ʾ�������Ӷ����ٳ��ȣ�����ֻҪlrc_n��20����Լ��Ҫ8M�ڴ棩���ڣ�һ���ǿ��Խ��ܵĴ�С
 *encode_offset �洢ÿһ�����ݿ��ڱ��������Ӧ��encode_task����opeλ�á� Ҳ����encode_offset[i] �� encode_offset[i+1] -1 ��Ӧ��ope���õ����ݿ�i
 */
static task *encode_task;
static task **decode_task_cache;
static int encode_offset[lrc_k + 1];

/*
 *temp_var ��Ӧ�����ʼ����Ҫ��һЩ������lrc_init()����free_temp_var()�ͻ��ͷ���Щ�ڴ档
 *mult_table �洢GF(2^8)�ĳ˷���
 *div_table �洢GF(2^8)�ĳ�����
 *encode_matrix �洢��ʼ�ı������
 *bitmap[k]�洢һ��Ϊ��ά���󣬱�ʾk�������GF(2^8)�ж�Ӧ��bitmap
 *tot_bit[k]��ʾk������ݶ�Ӧ��bitmap���ж��ٸ�1������ͳ�Ƹ���������ռ䣩
 *block_task[k]��ʾһ��k���ݿ����ʱ��task�� ����encode_task ��������block_task�ĺϲ�
 *now_ope��Ӧÿ�������е�bitλ�����Ϊ0����ʾ���ʱ��Ӧ�ø��ƣ����������xor����
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
 *����GF(2^8)�ĳ˷��ͳ�����
 *��Դ����ʽΪx^8 + x^4 + x^3 + x^2 + 1  ���������ʽֱ������ http://blog.csdn.net/mengboy/article/details/1514445 ��ƪ���ͣ�
 *a��Ϊ�˼��㷽���ʼ��a^i��Ӧ��ֵ
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
			//����ʽ�˷�
			for (k = 0; k < 8; k++)
				if (j & (1 << k))  z ^= (i << k);
			mult_table[i][j] = 0;
			//����ʽ�˷����ȡ��(��λ����������2**8�ĸ�λ����ȡ�ࣩ
			for (k = 0; k < 15; k++)
				mult_table[i][j] ^= ((z & (1 << k)) == 0) ? 0 : a[k];
			//��⽻���ɽ��
			if (i > j && mult_table[i][j] != mult_table[j][i]) {
				fprintf(stderr, "error: %d %d %d %d\n", i, j, mult_table[i][j], mult_table[j][i]);
				exit(-1);
			}

			//�ӳ˷����еõ�������
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
 *��now_ope�����������Ϊ0��
 */
static void reset_now_ope(temp_var *var) {
	int i, j;
	for (i = 0; i < lrc_n; i++)
		for (j = 0; j < 8; j++)
			var->now_ope[i][j] = 0;
}


/*
 *��data[d] = data[d] xor (k * data[s]) ��������bitmap�µ��������У����洢��now����
 *��bitmap������� data[d][i] = sum(foreach(j) (k[i][j]*data[s][j])), GF(2^8)���棬sum ������ͬ��xor����
 *��Ϊk[i][j] ֻ��0��1����������ҵ�Ϊ0��0*data[s][j] == 0  ==> data[d][i]���䣬����0��������Ժ���
 *��k[i][j]Ϊ1 ��ʱ�� k[i][j]*data[s][j] == data[s][j] �����Ͳ���Ҫ�˷���ֱ�Ӽ���sum, Ҳ����xor��
 *����Ϊһ��ʼdata[d][i]��Ҫ��ֵΪ0����ô���൱��һ������k[i][j] == 1ʱ��data[d][i] = data[s][j]; ����opeΪ0�����˺���������Ӧ��xor����
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
 *���ȹ������MS ������ķ���������������
 *Ȼ��Ѿ���ת��Ϊһ��������
 *ͬʱ������block_task
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
	//����ope��Դ��ַ�������򣬲�������Ӧ��encode_offset, �Ա��ڶ���һ��������src block���ж�β���������CPU cache locality
	//����tot_ope�����󡣾�ֱ�ӽ���ѡ������
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
 *�������档�����ŵ�inv����
 *�㷨��
 *    ��Ҫ����ľ���A �������� AE ��EΪ��λ����Ȼ������м�����㣬ת��Ϊ EB ������ô����B��ΪA����
 */
static int invert_matrix(temp_var *var, int *mat, int *inv, int rows) {
	int i, j, k, x, rs2, rs, temp, inverse;
	gf **mult_table = var->mult_table;
	gf **div_table = var->div_table;
	k = 0;
	//��ʼ��ΪE��
	for (i = 0; i < rows; i++) {
		for (j = 0; j < rows; j++) {
			inv[k] = (i == j) ? 1 : 0;
			k++;
		}
	}

	//��ÿһ�н��в���
	for (i = 0; i < rows; i++) {
		rs = rows * i;

		//���mat[i][i] Ϊ0���Ǿ��ҳ�mat[j][i]��Ϊ0��j�����ӵ�i�С� ����Ҳ�������˵���þ��������档�˳�
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
		//��mat[i][i]��Ϊ1��
		temp = mat[k];
		if (temp != 1) {
			inverse = div(1, temp);
			for (j = 0; j < rows; j++) {
				mat[rs + j] = mult(mat[rs + j], inverse);
				inv[rs + j] = mult(inv[rs + j], inverse);
			}
		}

		//��ȥmat[j][i] (j!=i)��ֵ
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
 *�����ṩ�Ĳ�������������ɽ������
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
 *�����ṩ��ѡ���к��е���Ϣ�������½�����󣬲����ɶ�Ӧ�Ĳ�����
 *need_decode��Ҫ����ı����Ϣ����ʵ�����𻵿��state
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

	//���ݱ���������������
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
 *����state�������Ҫ�Ĳ��������ж��ܷ���롣
 *����ܽ���͹�����Ӧ��task
 *need_para ��Ҫ�����ݿ��ǡ�
 *selected_group ��ЩС�鱻ѡ��ı��
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

	// �����ҳ���Ҫ�Ĵ���������飬��ͳ����Ҫȫ�ֱ����ĸ���    
	k = 0;
	temp = 0;  // ͳ����Ҫ���ٸ�ȫ�ֱ����
	for (i = 0; i < lrc_l; i++) {
		int b = 0;  // b�൱��һ��bool���ͣ����������������������ݣ���ô��bΪtrue������Ϊfalse
		for (j = 0; j < lrc_group_size; j++) b |= state[k + j];
		b |= (state[lrc_k + i]);
		selected_group[i] = b;  // ��������ݿ��𻵣���ô��ʾ�������������Ҫ������ġ�
		if (b) {
			for (j = 0; j < lrc_group_size; j++)
				if (state[k + j] == 0) {      //����������û���𻵣���ô������ݾ���Ҫ���ṩ������
					need_para[k + j] = 1;
				}
				else {
					if (need_para[lrc_k + i] == 0 && state[lrc_k + i] == 0) {   //���С������û���𻵣�����û�б����Ϊ��Ҫ���ݿ飬��ô�Ϳ�����С������������˳�����ݿ�
						need_para[lrc_k + i] = 1;
					}
					else {   //������Ҫһ��ȫ�����ݿ�ſ��Խ���
						temp++;
					}
				}
		}
		k += lrc_group_size;
	}
	for (i = lrc_k + lrc_l; i < lrc_n; i++) {   //�����һ��ȫ�ֱ���鱻˳�����ǵ���Ҫһ��ȫ�����ݿ�
		if (state[i] == 1) temp++;
	}

	// �����Ҫ��ȫ�ֱ�����������ȫ�ֱ������������ܱ����롣ֱ���˳���
	if (temp > lrc_r) {
		free(need_para);
		free(selected_group);
		return 0;
	}

	// �����Ҫȫ�����ݿ���룬�Ǿ���ζ������С������ݿ鶼��Ҫ�ṩ�����롣
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

	//������Ҫ�Ŀռ�
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

	// ������Ҫ���к���
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

	temp = 0; //���ٸ��������
	for (i = lrc_k; i < lrc_n; i++) {
		if (state[i] == 1) {
			temp += block_task[i].tot_ope;
		}
	}

	//����б�����𻵣���ôֻҪֱ�Ӷ�ȡ֮ǰԤ�����block_task,�ϲ����ղż����task���棬�Ϳ�����ɽ��롣
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
 *�������ݿ��𻵵������ ʹ������
 *stateΪһ������lrc_n�����飬1��ʾ��Ӧ�����ݿ��Ѿ�����
 *seat��ʾĿǰ�������ڼ���
 *tot��ʾĿǰһ���ж��ٿ����𻵵ġ� tot���ܳ���lrc_l+lrc_r
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
 *����decode_task_cache
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
 *��ʼ��bitmap,tot_bit,������now_ope�Ŀռ�
 *kת���ɶ�Ӧ��bitmap������㷨�ο�Jerasure��
 *�����㷨�� bitmap[k][i] ��Ӧk*2^i�Ķ�����λ
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
 *�ͷ�var���ڴ档 ��Ҫ����һ��var��������ָ�룬��֤�����ͷ�
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
 *���ݸ�����operation���У������ݽ��б���롣
 *�������ٶȵĺ��ġ�
 *��packetΪ��λ���б���롣
 *����һ��packet������Ϊ8�飬��Ϊb0,b2,...,b7���ֱ��Ӧ������ݵ�bitλ�����ݡ�
 *(b0i,b1i,...,b7i)��Ӧpacket�ĵ�i������ ��b0i��ʾb0��i��bitλ����
 * ���ݾ���Ĳ����ԣ��Ϳ��Խ��м������㡣ÿ��ȡ��sizeof��long��������һ�����㡣��Ҳ����cpu�ֳ���
 * �������ݷ������и���packetΪ��λ��������ÿ��byte���У���˲�ͬ��packet_size�ᵼ�²�ͬ�ļ�������
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
	lrc_free();              // ����֮ǰ�ڴ�û���ͷţ�������ڴ�й¶��

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











