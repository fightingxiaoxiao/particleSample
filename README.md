# particleSample

## 简介

用于对OpenFOAM的lagrangian颗粒进行采样。

## 安装

```shell
cd $HOME/OpenFOAM/OpenFOAM-<version>
mkdir extend
cd extend
git clone https://github.com/fightingxiaoxiao/particleSample

cd $HOME/OpenFOAM/OpenFOAM-<version>/applications/utilities/postProcessing/lagrangian
cp -r $HOME/OpenFOAM/OpenFOAM-<version>/extend/particleSample .

cd particleSample
wmake

```

## 使用

将[particleSampleProperties](particleSampleProperties)放入算例的constant目录。

在算例目录执行`particleSample`。如果是并行算例，执行`mpirun -n <number> particleSample -parallel`。
