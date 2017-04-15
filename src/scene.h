#pragma once
#ifndef _SCENE_H_
#define _SCENE_H_

#include"medium.h"
#include"material.h"
#include"common.h"


// 媒質のリスト
//Medium milkMedium(Color(4.5513, 5.8294, 7.136), Color(0.0015333, 0.0046, 0.019933), 1.3);
Medium milkMedium(Color(0.45513, 0.58294, 0.7136), Color(0.00015333, 0.00046, 0.0019933), 1.3);
Medium cokeMedium(Color(8.9053e-05, 8.372e-05, 0) * 1.0, Color(0.10014, 0.16503, 0.2468) * 1.0, 1.3);
Medium diamondMedium(Color(EPS, EPS, EPS) * 1.0, Color(EPS, EPS, EPS) * 1.0, 1.5);
// マテリアルのリスト
Material lightMat(Color(36.0), Color(0.0), DIFFUSE);//TODO弱めた
Material grayMat(Color(0.0), Color(0.75), DIFFUSE);
Material redMat(Color(0.0), Color(0.75, 0.25, 0.25), DIFFUSE);
Material greenMat(Color(0.0), Color(0.25, 0.75, 0.25), DIFFUSE);
Material blueMat(Color(0.0), Color(0.25, 0.25, 0.75), DIFFUSE);
Material mirrorMat(Color(0.0), Color(0.999), SPECULAR);
Material glassMat(Color(0.0), Color(0.999), REFRACTION);
Material milkMat(Color(0.0), Color(0.999), TRANSLUCENT, milkMedium);
Material cokeMat(Color(0.0), Color(0.999), TRANSLUCENT, cokeMedium);
Material diamondMat(Color(0.0), Color(0.999), REFRACTION, diamondMedium);




#endif _SCENE_H_