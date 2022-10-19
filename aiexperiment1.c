#include "stdio.h"
#include "stdlib.h"
#include "time.h"
#define cityNum 10
#define popSize 10
#define croRate 0.85
#define mutRate 0.1
#define MAX 999

//定义染色体的结构
struct Chrom
{
    int cityArr[cityNum]; //向量一共十个城市应该表示的是顺序
    char name;            //染色体名称
    float adapt;
    int dis; //距离，也就是路径长
};
struct Chrom genes[popSize];    //染色体组，数组内容十个
struct Chrom genesNew[popSize]; //用来存放新的染色体组
struct Chrom temp;              //临时的染色体数据

char names[cityNum] = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J'}; //十个城市的名称

int distance[cityNum][cityNum] = {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9},
                                  {1, 0, 1, 2, 3, 4, 5, 6, 7, 8},
                                  {2, 1, 0, 1, 2, 3, 4, 5, 6, 7},
                                  {3, 2, 1, 0, 1, 2, 3, 4, 5, 6},
                                  {4, 3, 2, 1, 0, 1, 2, 3, 4, 5},
                                  {5, 4, 3, 2, 1, 0, 1, 2, 3, 4},
                                  {6, 5, 4, 3, 2, 1, 0, 1, 2, 3},
                                  {7, 6, 5, 4, 3, 2, 1, 0, 1, 2},
                                  {8, 7, 6, 5, 4, 3, 2, 1, 0, 1},
                                  {9, 8, 7, 6, 5, 4, 3, 2, 1, 0}}; //城市之间的对应距离

void initGroup() //初始化随机数据获取最初的染色体组
{
    int i, j, k;
    int t = 0;
    int flag = 0;
    srand(time(NULL));            //初始化随机数生成器
    for (i = 0; i < popSize; i++) //最外层循环表示的是genes[i]也就是一开始的染色体组数据
    {

        temp.name = names[i]; //设置染色体的名称
        temp.adapt = 0.0f;    //初始化适应值也是健康值
        temp.dis = 0;

        for (j = 0; j < cityNum;) //把十个城市依次放入获取genes[i].cityArr[cityNum]
        {
            t = rand() % cityNum; //随机生成一个向量当作是行走路线
            flag = 1;             //把标识位置设置为1
            for (k = 0; k < j; k++)
            {
                if (genes[i].cityArr[k] == t) //不能重复必须是不重复的城市
                {
                    flag = 0;
                    break;
                }
            }
            if (flag) //如果是不重复的就放入
            {
                temp.cityArr[j] = t;
                genes[i] = temp;
                j++;
            }
        }
    }
}

void popFitness() //健康程度测量
{
    int i, n1, n2;
    for (i = 0; i < popSize; i++) //用来遍历染色体组
    {
        genes[i].dis = 0;
        for (int j = 1; j < cityNum; j++) //获取对应的总路线长度
        {
            n1 = genes[i].cityArr[j - 1];
            n2 = genes[i].cityArr[j];
            genes[i].dis += distance[n1][n2];
        }
        genes[i].dis += distance[genes[i].cityArr[0]][genes[i].cityArr[cityNum - 1]]; //再加上首末城市距离
        genes[i].adapt = (float)1 / genes[i].dis;                                     //得到健康度
    }
}

int chooseBest()
{
    int choose = 0;
    float best = 0.0f;
    best = genes[0].adapt; //从第0个开始
    for (int i = 0; i < popSize; i++)
    {
        if (genes[i].adapt < best)
        {
            best = genes[i].adapt;
            choose = i;
        }
    }
    return choose; //返回当前染色体组中最健康的染色体
}

void select()
{
    float biggestSum = 0.0f; //初始化一个健康度的综合
    float adapt_pro[popSize];
    float pick = 0.0f;
    int i;
    for (i = 0; i < popSize; i++) //遍历当前染色体组，获取健康度的总和
    {
        biggestSum += genes[i].adapt;
    }
    float a=0;
    for (i = 0; i < popSize; i++)
    {
        a=a+genes[i].adapt / biggestSum;
        adapt_pro[i] = a //获取概率占比并创建轮盘
    }
    for (i = 0; i < popSize; i++)
    {
        pick = (float)rand() / RAND_MAX; //所谓赌轮盘，也就是说随机出来的
                                         
        /********** Begin **********/

        for (int k = 0; k < popSize; k++)
        {
            a+=adapt_pro[k];
            if(pick<adapt_pro[k]){
                genesNew[i]=genes[k];
                break;
            }
        }

        /********** End **********/
    }
    for (i = 0; i < popSize; i++)
    {
        genes[i] = genesNew[i]; //更新染色体组
    }
}

void cross()
{
    float pick;
    int choice1, choice2;
    int pos1, pos2;
    int temp;
    int index1, index2;
    int move = 0;
    while (move < popSize - 1)
    {
        pick = (float)rand() / RAND_MAX;
        if (pick > croRate) //如果比交叉概率大就把move移动两位，因为choice1和choice2相邻
        {
            move += 2;
            continue;
        }
        choice1 = move;
        choice2 = move + 1;
        pos1 = rand() % popSize; //随机生成pos1和pos2
        pos2 = rand() % popSize;
        while (pos1 > popSize - 2 || pos1 < 1) //生成不对要重新生成
        {
            pos1 = rand() % popSize;
        }
        while (pos2 > popSize - 2 || pos2 < 1)
        {
            pos2 = rand() % popSize;
        }

        if (pos1 > pos2) //要保证pos2大于pos1
        {
            temp = pos1;
            pos1 = pos2;
            pos2 = temp;
        }

        for (int j = pos1; j <= pos2; j++) //交换相邻染色体的部分片段
        {
            temp = genes[choice1].cityArr[j];
            genes[choice1].cityArr[j] = genes[choice2].cityArr[j];
            genes[choice2].cityArr[j] = temp;
        }

        num1 = 0;
        num2 = 0;

        if (pos1 > 0 && pos2 < popSize - 1) //必须检查是否冲突
        {
            /********** Begin **********/
            int flag1=0;
            while(1){
                int flag2=0;
                for (int i=pos1;i<=pos2;i++){
                    if(genes[choice1].cityArr[flag1]==genes[choice1].cityArr[i]){
                        genes[choice1].cityArr[flag1]=genes[choice2].cityArr[i];
                        flag2++;
                    }
                    if(genes[choice2].cityArr[flag1]==genes[choice2].cityArr[i]){
                        genes[choice2].cityArr[flag1]=genes[choice1].cityArr[i];
                        flag2++;
                    }
                }
                if(flag2==0){
                    flag1++;
                }
                if(flag1==pos1)break;
            }

            /********** End **********/

            flag1=pos2+1;
            while (1)
            {
                int flag2=0;
                for (int i=pos1;i<=pos2;i++){
                    if(genes[choice1].cityArr[flag1]==genes[choice1].cityArr[i]){
                        genes[choice1].cityArr[flag1]=genes[choice2].cityArr[i];
                        flag2++;
                    }
                    if(genes[choice2].cityArr[flag1]==genes[choice2].cityArr[i]){
                        genes[choice2].cityArr[flag1]=genes[choice1].cityArr[i];
                        flag2++;
                    }
                }
                if(flag2==0){
                    flag1++;
                }
                if(flag1==cityNum)break;
            }
            
        }

        move += 2;
    }
}

void mutation() //突变函数
{
    double pick;
    int pos1, pos2, temp;
    for (int i = 0; i < popSize; i++)
    {
        pick = (float)rand() / RAND_MAX;
        if (pick > mutRate)
        {
            continue;
        }
        pos1 = rand() % popSize;
        pos2 = rand() % popSize;
        while (pos1 > popSize - 1)
        {
            pos1 = rand() % popSize;
        }
        while (pos2 > popSize - 1)
        {
            pos2 = rand() % popSize;
        }

        int a = genes[i].dis;
        temp = genes[i].cityArr[pos1];
        genes[i].cityArr[pos1] = genes[i].cityArr[pos2];
        genes[i].cityArr[pos2] = temp;

        popFitness();
        if (genes[i].dis > a)
        {
            temp = genes[i].cityArr[pos1];
            genes[i].cityArr[pos1] = genes[i].cityArr[pos2];
            genes[i].cityArr[pos2] = temp;
        }
    }
}
