/*
 *实现论文 “Erasure Coding in Windows Azure Storage” 里面介绍的Local Reconstruction Codes算法
 *Author：黄曙权（shughuangqqq@gmail.com)
 */


#ifndef _LRC_H
#define _LRC_H


 /*
  *		定义Local Reconstruction Codes 算法的基本参数
  *		lrc_k 表示数据块个数
  *		lrc_l 表示分成的小组个数，由于采用GF（2^8) 只能为 2  若为3，那么每组长度只能为3.
  *		lrc_r 表示全局编码块个数
  *		lrc_n 所有块个数 也就是 lrc_k+lrc_l+lrc_n
  *		lrc_group_size 表示每个小组有多少个数据块，也就是 lrc_k / lrc_l;
  */

#define lrc_k 12      
#define lrc_l 2              
#define lrc_r 2                
#define lrc_n 16
#define lrc_group_size 6

  /*
   *		lrc_packet_size 表示处理数据的基本大小单位，
   *						注意，如果编码和解码时候的packet_size不同，会导致解码无效。
   *						同时由于cpu优化等，packet_size不能太小，推荐在1024以上。
   *		lrc_bitpacket_size 表示在一个基本单位中每个bit位长度。
   *							由于取GF(2^8),所以每个bit的长度为packet_size / 8。
   *							注意，由于采用cpu字长进行优化，bitpacket_size  必须是sizeof(long)的倍数。（也就是计算机字长）
   *		lrc_bitpacket_skip  log(bitpacket_size) 为了程序优化而设置。
   */
#define lrc_packet_size 4096   
#define lrc_bitpacket_size 512    
#define lrc_bitpacket_skip 9      


   /*
	* 错误信息
	*/
#define LRC_ERROR_OUT_OF_MEMORY (-2)


typedef unsigned char gf;


/*
 *lrc初始化，生成编码和解码需要的操作串。
 *可能因为内存申请失败而退出。
 *正确初始化返回0；
 */
int lrc_init();

/*
 *释放lrc的内存。清除由lrc_init()产生的编解码信息。
 */
void lrc_free();

/*
 * 传入数据data。长度为lrc_n的数组指针。依次表示数据块，小组编码块，全局编码块的数据。
 * 再解码时，data前lrc_k为给定数据，编码后会把编码快一次放入 后（lrc_l+lrc_r)块。
 * size 为每个数据块长度
 * 注意：后（lrc_l + lrc_r)块的空间需要实现申请，本方法不处理空间申请和销毁操作。
 *		size必须是lrc_packet_size的倍数
 */
void lrc_encoding(gf **data, long size);

/*
 *只有部分数据更新的时候，更新编码块
 *newData 为长度为lrc_n的数组，依次表示数据块，小组编码块，全局编码块。如果数据块不为NULL，表示对应块有数据更新，并且为编码块申请好空间,产生的编码放入对应位置。
 *oldData 把旧数据和旧的编码块存放在对应数据，位置和newData的位置一一对应，这个方法不会修改oldData
 */
void lrc_encoding_part(gf **newData, gf **oldData, long size);


/*
 *提供一个数据块损坏状态，返回一个标示符，用于参数请求和解码使用。
 *state 为长度为lrc_n的数组，依次表示数据块，小组编码块，全局编码块。如果为1，表示块损坏，需要解码
 */
long lrc_state_to_long(int *state);

/*
 *返回需要数据状态。
 *返回值为长度为lrc_n的数组，依次表示数据块，小组编码块，全局编码块。如果为1，解码时需要提供的数据。
 *如果不能编码，返回NULL
 *注意，本返回数组为只读。不能修改，不能销毁。
 */
int *lrc_need_para(long type);


/*
 *解码。
 *提供必要的data。
 *长度为lrc_n的数组指针。依次表示数据块，小组编码块，全局编码块的数据。
 *需要的数据放入对应的位置，并为需要解码的块申请空间。而不需要的块可以赋值为NULL。
 *size 为每个数据块长度，且必须是lrc_packet_size的倍数
 */
void lrc_decoding(long type, gf **data, long size);


/*
 *注意：
 *		lrc_encoding,lrc_state_to_long,lrc_need_para,lrc_decoding 为线程安全，可以多线程同时调用。
 *		为了加快数据编码，所有的方法不会对传入数据合法性做检查，需要自己保证参数合法，不然会引起崩溃。
 */

#endif
