/*
 *ʵ������ ��Erasure Coding in Windows Azure Storage�� ������ܵ�Local Reconstruction Codes�㷨
 *Author������Ȩ��shughuangqqq@gmail.com)
 */


#ifndef _LRC_H
#define _LRC_H


 /*
  *		����Local Reconstruction Codes �㷨�Ļ�������
  *		lrc_k ��ʾ���ݿ����
  *		lrc_l ��ʾ�ֳɵ�С����������ڲ���GF��2^8) ֻ��Ϊ 2  ��Ϊ3����ôÿ�鳤��ֻ��Ϊ3.
  *		lrc_r ��ʾȫ�ֱ�������
  *		lrc_n ���п���� Ҳ���� lrc_k+lrc_l+lrc_n
  *		lrc_group_size ��ʾÿ��С���ж��ٸ����ݿ飬Ҳ���� lrc_k / lrc_l;
  */

#define lrc_k 12      
#define lrc_l 2              
#define lrc_r 2                
#define lrc_n 16
#define lrc_group_size 6

  /*
   *		lrc_packet_size ��ʾ�������ݵĻ�����С��λ��
   *						ע�⣬�������ͽ���ʱ���packet_size��ͬ���ᵼ�½�����Ч��
   *						ͬʱ����cpu�Ż��ȣ�packet_size����̫С���Ƽ���1024���ϡ�
   *		lrc_bitpacket_size ��ʾ��һ��������λ��ÿ��bitλ���ȡ�
   *							����ȡGF(2^8),����ÿ��bit�ĳ���Ϊpacket_size / 8��
   *							ע�⣬���ڲ���cpu�ֳ������Ż���bitpacket_size  ������sizeof(long)�ı�������Ҳ���Ǽ�����ֳ���
   *		lrc_bitpacket_skip  log(bitpacket_size) Ϊ�˳����Ż������á�
   */
#define lrc_packet_size 4096   
#define lrc_bitpacket_size 512    
#define lrc_bitpacket_skip 9      


   /*
	* ������Ϣ
	*/
#define LRC_ERROR_OUT_OF_MEMORY (-2)


typedef unsigned char gf;


/*
 *lrc��ʼ�������ɱ���ͽ�����Ҫ�Ĳ�������
 *������Ϊ�ڴ�����ʧ�ܶ��˳���
 *��ȷ��ʼ������0��
 */
int lrc_init();

/*
 *�ͷ�lrc���ڴ档�����lrc_init()�����ı������Ϣ��
 */
void lrc_free();

/*
 * ��������data������Ϊlrc_n������ָ�롣���α�ʾ���ݿ飬С�����飬ȫ�ֱ��������ݡ�
 * �ٽ���ʱ��dataǰlrc_kΪ�������ݣ�������ѱ����һ�η��� ��lrc_l+lrc_r)�顣
 * size Ϊÿ�����ݿ鳤��
 * ע�⣺��lrc_l + lrc_r)��Ŀռ���Ҫʵ�����룬������������ռ���������ٲ�����
 *		size������lrc_packet_size�ı���
 */
void lrc_encoding(gf **data, long size);

/*
 *ֻ�в������ݸ��µ�ʱ�򣬸��±����
 *newData Ϊ����Ϊlrc_n�����飬���α�ʾ���ݿ飬С�����飬ȫ�ֱ���顣������ݿ鲻ΪNULL����ʾ��Ӧ�������ݸ��£�����Ϊ���������ÿռ�,�����ı�������Ӧλ�á�
 *oldData �Ѿ����ݺ;ɵı�������ڶ�Ӧ���ݣ�λ�ú�newData��λ��һһ��Ӧ��������������޸�oldData
 */
void lrc_encoding_part(gf **newData, gf **oldData, long size);


/*
 *�ṩһ�����ݿ���״̬������һ����ʾ�������ڲ�������ͽ���ʹ�á�
 *state Ϊ����Ϊlrc_n�����飬���α�ʾ���ݿ飬С�����飬ȫ�ֱ���顣���Ϊ1����ʾ���𻵣���Ҫ����
 */
long lrc_state_to_long(int *state);

/*
 *������Ҫ����״̬��
 *����ֵΪ����Ϊlrc_n�����飬���α�ʾ���ݿ飬С�����飬ȫ�ֱ���顣���Ϊ1������ʱ��Ҫ�ṩ�����ݡ�
 *������ܱ��룬����NULL
 *ע�⣬����������Ϊֻ���������޸ģ��������١�
 */
int *lrc_need_para(long type);


/*
 *���롣
 *�ṩ��Ҫ��data��
 *����Ϊlrc_n������ָ�롣���α�ʾ���ݿ飬С�����飬ȫ�ֱ��������ݡ�
 *��Ҫ�����ݷ����Ӧ��λ�ã���Ϊ��Ҫ����Ŀ�����ռ䡣������Ҫ�Ŀ���Ը�ֵΪNULL��
 *size Ϊÿ�����ݿ鳤�ȣ��ұ�����lrc_packet_size�ı���
 */
void lrc_decoding(long type, gf **data, long size);


/*
 *ע�⣺
 *		lrc_encoding,lrc_state_to_long,lrc_need_para,lrc_decoding Ϊ�̰߳�ȫ�����Զ��߳�ͬʱ���á�
 *		Ϊ�˼ӿ����ݱ��룬���еķ�������Դ������ݺϷ�������飬��Ҫ�Լ���֤�����Ϸ�����Ȼ�����������
 */

#endif
