# easyRATE

Eyring equation 속의 $\kappa$는 transmission coefficient로 터널링 효과를 고려한 상수입니다. 하지만 많은 경우 이 값을 알기 어렵기 때문에 1로 가정합니다. Transmission coefficient 값을 추정하는 방법은 크게 3가지 정도 알려져 있는데 Wigner, Skodje, Eckart 선생님들이 각각 만드신 방법들이 있습니다. 이중 Wigner 방법을 제외한 두 방법은 비교적 많은 정보가 필요합니다. easyRATE는 Wigner transmission correction 방법을 이용해 최소한의 정보로 간단한 transmission correction을 수행해볼 수 있는 코드입니다.

Wigner transmission correction은 오직 TS에서의 `Imaginary Vibrational Frequency`에 대한 정보만을 필요로 합니다. 이 외에 Rate constant를 구하기 위해서 `Reactant와 TS에서의 Gibbs Energy`, 그리고 `온도`를 필요로 합니다. 이외에 옵션으로 Reaction Symmetry 보정을 위한 `Reactant와 TS의 xyz file`이 필요합니다.

* 필요 정보
  * 온도
  * 전이상태에서의 허수 진동수
  * 깁스 에너지
  * xyz 좌표 파일 (option)

## How to use
두가지 사용 모드가 존재합니다.

명령행 인수가 있으면 input mode로 실행되고 없으면 interactive mode로 실행됩니다.

1. input mode  :  input 파일을 작성한 후 커맨드줄 인수로 제공
2. interactive mode  :  커맨드 라인에서 필요한 값들을 직접 제공함

```
# interactive mode #

> python easyRATE.py
```
```
# input mode #

> python easyRATE.py input.inp
```
필요한 인풋 파일의 구조는 다음과 같습니다. ( 작성 방법은 input.inp 파일에 더 자세히 적어뒀습니다. )
```
# Imaginary Frequency in TS #
Imaginary Frequency [cm^-1] = -1530


# Temperatures #
Temperatures list [K] = 300 500 10


# Special Temperatures (optional) #
Special Temperatures [K] = auto


# Unit Type #
#   1   :  Hartree  
#   2   :  kcal/mol
#   3   :  kJ/mol
Unit Type = 2


# Thermodynamic Energy #
Gibbs(Reactant) = -74
Gibbs(TS) = -39


# RXN symmetry number #
RXN Symmetry Number = 1


# consumption time # 
consumption_percent = 50 


# XYZ file (optional) #
Reactant xyz file = 
TS xyz file = 
```


## Equations

* wigner transmission correction
```math
\kappa_{wig}(T) = 1 + \frac{1}{24}{\left( \frac{h |v^{\ddagger}|}{k_{b} T} \right)} ^2
```
* Eyring equation
```math
k(T) = \kappa_{wig}(T) \sigma_{sym} \frac{k_b T}{h} e^{-\Delta G^{\ddagger} /RT }
```
* Reaction Symmetry Number
```math
\sigma_{sym} = \frac{\sigma_{rot, Rc}}{\sigma_{rot, TS}}
```

## Limitation

* Wigner correction의 경우 Eckart correction과 비교했을 때 저온에서 부정확함
* gas phase rxn에 대해서만 지원함

## Todo

- [ ] Reaction symmetry auto-estimation에 필요한 point 그룹별 roatational symmetry number 추가
- [ ] 현재 A -> TS -> B 반응만 지원함. 더 많은 반응 지원

## Reference
* Int J Quantum Chem. (2018) ; e25686.
* Theor Chem Account (2007) 118, 813–826
