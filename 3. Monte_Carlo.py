#https://zh.wikipedia.org/wiki/%E8%92%99%E5%9C%B0%E5%8D%A1%E7%BE%85%E6%96%B9%E6%B3%95
#https://zh.wikipedia.org/wiki/%E8%92%99%E7%89%B9%E5%8D%A1%E6%B4%9B%E6%A0%91%E6%90%9C%E7%B4%A2
#在面積為1的正方形中有1/4正圓，亂數選擇點位置，點在兩者面積比值
#1:pi/4 = 全部點數量:圓內點數量  =>  pi = 圓內/全部 *4

import random as rd
def cal_pi():
    count = 0
    total = 1000000
    for i in range(total):
        x = rd.random()
        y = rd.random()
        if (x**2 + y**2 <= 1):
            count +=1
    return format(count/total *4, ".20f")
print(cal_pi())

