# easyRATE

## How to use
To be add
```
python easyRATE.py
```

## Example output
```
Temperature       kappa_wig             k_eyring               dacay time
----------------------------------------------------------------------------------------------
273.15            3.6011E+00            2.0320E-19            108164677506.0 years 9 months
293.15            3.2583E+00            3.0103E-17            730138261.0 years 5 months
298.15            3.1832E+00            9.4606E-17            232327028.0 years 11 months
400.00            2.2129E+00            2.5782E-09            8.0 years 6 months
410.00            2.1545E+00            8.7790E-09            2.0 years 6 months
420.00            2.1002E+00            2.8215E-08            9 months 14 days
430.00            2.0496E+00            8.5932E-08            3 months 3 days
440.00            2.0024E+00            2.4892E-07            1 months 2 days
450.00            1.9584E+00            6.8813E-07            11 days 15 hours
460.00            1.9172E+00            1.8209E-06            4 days 9 hours
470.00            1.8785E+00            4.6252E-06            1 days 17 hours
480.00            1.8423E+00            1.1306E-05            17 hours 1 minutes
490.00            1.8083E+00            2.6659E-05            7 hours 13 minutes
500.00            1.7763E+00            6.0767E-05            3 hours 10 minutes
```

## Equations

to be add

```math

\kappa_{wig}(T) = 1 + \frac{1}{24}{\left( \frac{h |v^{\ddagger}|}{k_{b} T} \right)} ^2

```

## Limitation

* Wigner correction의 경우 Eckart correction과 비교했을 때 저온에서 부정확함

## Todo

- [ ] Reaction symmetry auto-estimation code 추가
- [ ] Rot. symmetry number dictionary에 값 추가
- [ ] Input mode 추가
- [ ] output 파일 생성 코드 추가
- [ ] 현재 A -> TS -> B 반응만 지원함. 더 많은 반응 지원
- [ ] pyTUN mode 추가

## Reference
* Int J Quantum Chem. (2018) ; e25686.
* Theor Chem Account (2007) 118, 813–826
