# <font face="黑体">数值积分方法</font>

<font face="宋体">代码编写：御河DE天街</font><br/>
<font face="Times New Roman">Powered by: FantasySilence<br/>
<br/>
<font face="宋体">注意：对于瑕积分而言数值积分方法或部分不适用，可以在传入参数时，以分段函数的形式指定瑕点，从而增加适应性。</font><br/>

## <font face="黑体">1.自适应积分</font>

### <font face="黑体">1.1.初始化对象</font>
<font face="宋体">首先导入对应的模块<br/>

```
from Numerical_Intergration.AdaptiveInteg import AdaptiveIntergralAlgorithm
```

#### <font face="黑体">1.1.1.参数说明</font>

```
AdaptiveIntergralAlgorithm(int_fun, int_internal, eps=1e-8)
```
<font face="宋体">其中，```int_fun```为被积函数；```int_internal```为积分区间；```eps```为精度，默认为1e-8。

#### <font face="黑体">1.1.2.类方法(属性)简介</font>

<font face="宋体">计算积分：```self.cal_int()```<br/>
积分值：```self.int_value```<br/>
划分区间：```self.x_node```<br/>

### <font face="黑体">1.2.调用示例</font>

```
import numpy as np
from Numerical_Intergration.AdaptiveInteg import AdaptiveIntergralAlgorithm

def fun(x):
    return 1/(np.sin(6*np.pi*x)+np.sqrt(x))

ada = AdaptiveIntergralAlgorithm(fun, [1,2], eps=1e-5)
ada.cal_int()
print("自适应积分值：", ada.int_value)
print("划分区间：", ada.x_node)
print("划分节点数：", len(ada.x_node))
```

## <font face="黑体">2.复合求积公式</font>

### <font face="黑体">2.1.初始化对象</font>
<font face="宋体">首先导入对应的模块<br/>

```
from Numerical_Intergration.CompositeFormula import CompositQuadratureIntegration
```

#### <font face="黑体">2.1.1.参数说明</font>

```
CompositQuadratureIntegration(int_fun, int_internal,internal_num=16, int_type="simpson")
```

<font face="宋体">其中，```int_fun```为被积函数，```int_internal```为积分区间，```internal_num```为区间划分数量，默认为16，```int_type```为积分公式类型，默认为辛普森复合求积公式，此外还有复合梯形公式(```trapezoid```)和复合科特斯公式(```cotes```)。

#### <font face="黑体">2.1.2.类方法(属性)简介</font>

<font face="宋体">计算积分：```self.cal_int()```<br/>
积分值：```self.int_value```<br/>
积分余项：```self.int_remainder```<br/>

### <font face="黑体">2.2.调用示例</font>

```
import sympy as sp
from Numerical_Intergration.CompositeFormula import CompositQuadratureIntegratio

t =sp.Symbol('t')
int_fun = sp.sin(t)

print("="*50)
cqi1 = CompositQuadratureIntegration(int_fun, [0, np.pi/2], int_type="trapezoid")
int_value1 = cqi1.cal_int()
print("复合梯形公式积分值：", int_value1)
print("复合梯形公式积分的余项：",cqi1.int_remainder)
print("="*50)
cqi2 = CompositQuadratureIntegration(int_fun, [0, np.pi/2], int_type="simpson")
int_value2 = cqi2.cal_int()
print("复合辛普森公式积分值：", int_value2)
print("复合辛普森公式积分的余项：",cqi2.int_remainder)
print("="*50)
cqi3 = CompositQuadratureIntegration(int_fun, [0, np.pi/2], int_type="cotes")
int_value3 = cqi3.cal_int()
print("复合科特斯公式积分值：", int_value3)
print("复合科特斯公式积分的余项：",cqi3.int_remainder)
print("="*50)
```

## <font face="黑体">3.离散数据的B样条插值法积分</font>

### <font face="黑体">3.1.初始化对象</font>
<font face="宋体">首先导入对应的模块<br/>

```
from Numerical_Intergration.CubicBsplineIterpInt import CubicBsplineInterpolationIntergration
```

#### <font face="黑体">3.1.1.参数说明</font>

```
CubicBsplineInterpolationIntergration(x, y)
```

<font face="宋体">其中，```x```和```y```是离散数据点，离散数据点的长度需一致并大于3，并且数据点之间是等距的。

#### <font face="黑体">3.1.2.类方法(属性)简介</font>

<font face="宋体">计算积分：```self.cal_int()```<br/>
积分值：```self.int_value```<br/>

### <font face="黑体">3.2.调用示例</font>

```
import numpy as np
from Numerical_Intergration.CubicBsplineIterpInt import CubicBsplineInterpolationIntergration

x = np.linspace(0, 24, 25)
y = np.array([0, 0.45, 1.79, 4.02, 7.15, 11.18, 16.09, 21.90, 29.05,
                29.05, 29.05, 29.05, 29.05, 22.42, 17.9, 17.9, 17.9,
                17.9, 14.34, 11.01, 8.9, 6.54, 2.03, 0.55, 0])
cbsi1 = CubicBsplineInterpolationIntergration(x, y)
int_value1 = cbsi1.cal_int()
print(int_value1)
```

## <font face="黑体">4.高斯——切比雪夫求积公式</font>
### <font face="黑体">4.1.初始化对象</font>
<font face="宋体">首先导入对应的模块<br/>

```
from Numerical_Intergration.GuassChebyshevInt import GuassChebyshevIntergration
```

#### <font face="黑体">4.1.1.参数说明</font>

```
GuassChebyshevIntergration(int_fun, zeros_num=10, cheb_type=1)
```

<font face="宋体">其中，```int_fun```为被积函数，```zeros_num```为切比雪夫正交多项式零点个数，默认为10，```cheb_type```为切比雪夫类型，对应第一类或第二类切比雪夫多项式，默认为1(第一类切比雪夫多项式)

#### <font face="黑体">4.1.2.类方法(属性)简介</font>

<font face="宋体">计算积分：```self.cal_int()```<br/>
积分值：```self.int_value```<br/>
切比雪夫零点：```self.zero_points```<br/>
高斯型求积公式积分系数：```self.A_k```

### <font face="黑体">4.2.调用示例</font>

```
import numpy as np
from Numerical_Intergration.GuassChebyshevInt import GuassChebyshevIntergration

def fun1(x):
    # np.exp(x)/np.sqrt(1-x**2)
    return np.exp(x)/np.sqrt(1-x**2)

def fun2(x):
    # x**2*np.sqrt(1-x**2)
    return x**2*np.sqrt(1-x**2)

cheb1 = GuassChebyshevIntergration(fun1, zeros_num=10,cheb_type=1)
cheb1.cal_int()
print(cheb1.zero_points)
print(cheb1.A_k)
print("-"*60)
cheb2 = GuassChebyshevIntergration(fun2, zeros_num=10,cheb_type=2)
cheb2.cal_int()
print(cheb2.zero_points)
print(cheb2.A_k)
```

## <font face="黑体">5.高斯——埃尔米特求积公式</font>
### <font face="黑体">5.1.初始化对象</font>
<font face="宋体">首先导入对应的模块<br/>

```
from Numerical_Intergration.GuassHermiteInt import GuassHermiteIntergration
```

#### <font face="黑体">5.1.1.参数说明</font>

```
GuassHermiteIntergration(int_fun, zeros_num=10)
```

<font face="宋体">其中，```int_fun```为被积函数，```zeros_num```为埃尔米特正交多项式零点个数，默认为10

#### <font face="黑体">5.1.2.类方法(属性)简介</font>

<font face="宋体">计算积分：```self.cal_int()```<br/>
积分值：```self.int_value```<br/>
埃尔米特零点：```self.zero_points```<br/>
高斯型求积公式积分系数：```self.A_k```

### <font face="黑体">5.2.调用示例</font>

```
import numpy as np
from Numerical_Intergration.GuassHermiteInt import GuassHermiteIntergration

def fun(x):
    return x**2*np.exp(-x**2)

herm = GuassHermiteIntergration(fun, zeros_num=15)
herm.cal_int()
print(herm.zero_points)
print(herm.A_k)
print(herm.int_value)
```

## <font face="黑体">6.高斯——拉盖尔积分公式</font>
### <font face="黑体">6.1.初始化对象</font>
<font face="宋体">首先导入对应的模块<br/>

```
from Numerical_Intergration.GuassLaguerreInt import GuassLaguerreIntergration
```

#### <font face="黑体">6.1.1.参数说明</font>

```
GuassLaguerreIntergration(self, int_fun, int_internal, zeros_num=10)
```

<font face="宋体">其中，```int_fun```为被积函数，```int_internal```为积分区间([a, +∞)),```zeros_num```为拉盖尔正交多项式零点个数，默认为10

#### <font face="黑体">6.1.2.类方法(属性)简介</font>

<font face="宋体">计算积分：```self.cal_int()```<br/>
积分值：```self.int_value```<br/>
拉盖尔零点：```self.zero_points```<br/>
高斯型求积公式积分系数：```self.A_k```

### <font face="黑体">6.2.调用示例</font>

```
import numpy as np
from Numerical_Intergration.GuassLaguerreInt import GuassLaguerreIntergration

def fun(x):
    return np.sin(x)*np.exp(-x)

lagr = GuassLaguerreIntergration(fun, [-2, np.infty], zeros_num=15)
lagr.cal_int()
print(lagr.zero_points)
print(lagr.A_k)
print(lagr.int_value)
```

## <font face="黑体">7.龙贝格加速算法</font>
### <font face="黑体">7.1.初始化对象</font>
<font face="宋体">首先导入对应的模块<br/>

```
from Numerical_Intergration.RombergAccelerationInt import RombergAccelerationQuadrature
```

#### <font face="黑体">7.1.1.参数说明</font>

```
RombergAccelerationQuadrature(self, int_fun, int_internal, acceleration_num=10)
```

<font face="宋体">其中，```int_fun```为被积函数，```int_internal```为积分区间,```acceleration_num```为外推次数，默认为10

#### <font face="黑体">7.1.2.类方法(属性)简介</font>

<font face="宋体">计算积分：```self.cal_int()```<br/>
积分值：```self.int_value```<br/>
龙贝格加速计算表：```self.romberg_table```<br/>

### <font face="黑体">7.2.调用示例</font>

```
from Numerical_Intergration.RombergAccelerationInt import RombergAccelerationQuadrature

def fun(x):
    return x**(3/2)

raq1 = RombergAccelerationQuadrature(fun, [0, 1], acceleration_num=15)
raq1.cal_int()
print(raq1.int_value)
print(raq1.romberg_table)
```

## <font face="黑体">8.复合辛普森二重积分公式</font>
### <font face="黑体">8.1.初始化对象</font>
<font face="宋体">首先导入对应的模块<br/>

```
from Numerical_Intergration.CompositeDoubleSimpsonInt import CompositeDoubleSimpsonIntergration
```

#### <font face="黑体">8.1.1.参数说明</font>

```
CompositeDoubleSimpsonIntergration(self, int_fun, x_span, y_span, eps=1e-6, max_split=100, increment=10)
```

<font face="宋体">其中，```int_fun```为被积函数，```x_span```为```x```的积分区间，```y_span```为```y```的积分区间，```eps```为积分精度，默认为1e-6，```max_split```为最大划分次数,默认为100，```increment```为每次递增区间数,默认从10开始。

#### <font face="黑体">8.1.2.类方法(属性)简介</font>

<font face="宋体">计算积分：```self.cal_2d_int()```<br/>
绘制自适应积分过程精度收敛曲线：```self.plt_precision(is_show=True)```<br/>
划分区间数：```self.internal_num```<br/>
积分值：```self.int_value```<br/>

### <font face="黑体">8.2.调用示例</font>

```
import numpy as np
from Numerical_Intergration.CompositeDoubleSimpsonInt import CompositeDoubleSimpsonIntergration

def fun(x,y):
    return np.exp(-x**2-y**2)

cdsi = CompositeDoubleSimpsonIntergration(fun, [0,1], [0,1], eps=1e-15)
cdsi.cal_2d_int()
print("划分区间数:%d, 积分近似值:%.15f"%(cdsi.internal_num, cdsi.int_value))
cdsi.plt_precision()
```

## <font face="黑体">9.高斯-勒让德积分公式</font>
### <font face="黑体">9.1.初始化对象</font>
<font face="宋体">首先导入对应的模块<br/>

```
from Numerical_Intergration.GuassLegendreInt import GuassLegendreIntergration
```

#### <font face="黑体">9.1.1.参数说明</font>

```
GuassLegendreIntergration(self, int_fun, int_internal, zeros_num=10)
```

<font face="宋体">其中，<font face="宋体">其中，```int_fun```为被积函数，```int_internal```为积分区间,```zeros_num```为勒让德正交多项式零点个数，默认为10。

#### <font face="黑体">9.1.2.类方法(属性)简介</font>

<font face="宋体">计算积分：```self.cal_int()```<br/>
积分值：```self.int_value```<br/>
勒让德零点：```self.zero_points```<br/>
高斯型求积公式积分系数：```self.A_k```

### <font face="黑体">9.2.调用示例</font>

```
import numpy as np
from Numerical_Intergration.GuassLegendreInt import GuassLegendreIntergration

def fun(x):
    return np.sin(x)*np.exp(-x)
    
leg = GuassLegendreIntergration(fun, [0,8], zeros_num=10)
int_value = leg.cal_int()
print("积分值:%.15f"%int_value)
print(leg.A_k)
```

## <font face="黑体">10.高斯-勒让德二重积分公式</font>
### <font face="黑体">10.1.初始化对象</font>
<font face="宋体">首先导入对应的模块<br/>

```
from Numerical_Intergration.GuassLegendreDuobleInt import GuassLegendreDuobleIntergration
```

#### <font face="黑体">10.1.1.参数说明</font>

```
GuassLegendreDuobleIntergration(self, int_fun, x_span, y_span, zeros_num=10)
```

<font face="宋体">其中，```int_fun```为被积函数，```x_span```为```x```的积分区间，```y_span```为```y```的积分区间，```zeros_num```为勒让德正交多项式零点个数，默认为10。

#### <font face="黑体">10.1.2.类方法(属性)简介</font>

<font face="宋体">计算积分：```self.cal_2d_int()```<br/>
积分值：```self.int_value```<br/>

### <font face="黑体">10.2.调用示例</font>

```
import numpy as np
from Numerical_Intergration.GuassLegendreDuobleInt import GuassLegendreDuobleIntergration

def fun(x,y):
    return np.exp(-x**2-y**2)

gldi = GuassLegendreDuobleIntergration(fun, [0,1], [0,1], zeros_num=15)
res = gldi.cal_2d_int()
print("二重积分结果：", res)
```

## <font face="黑体">11.高斯-勒让德三重积分公式</font>
### <font face="黑体">11.1.初始化对象</font>
<font face="宋体">首先导入对应的模块<br/>

```
from Numerical_Intergration.GuassLegendreTripleInt import GuassLegendreTripleIntergration
```

#### <font face="黑体">11.1.1.参数说明</font>

```
GuassLegendreTripleIntergration(self, int_fun, x_span, y_span, z_span, zeros_num=None)
```

<font face="宋体">其中，```int_fun```为被积函数，```x_span```为```x```的积分区间，```y_span```为```y```的积分区间，，```z_span```为```z```的积分区间```zeros_num```为勒让德正交多项式零点个数。

#### <font face="黑体">11.1.2.类方法(属性)简介</font>

<font face="宋体">计算积分：```self.cal_3d_int()```<br/>
积分值：```self.int_value```<br/>

### <font face="黑体">11.2.调用示例</font>

```
import numpy as np
from Numerical_Intergration.GuassLegendreTripleInt import GuassLegendreTripleIntergration

def fun(x,y,z):
    return 4*x*z*np.exp(-x**2*y-z**2)

glti = GuassLegendreTripleIntergration(fun, [0,1], [0,np.pi], [0,np.pi], zeros_num=[11, 12, 15])
res = glti.cal_3d_int()
print("积分值：%.15f"%res)
```
