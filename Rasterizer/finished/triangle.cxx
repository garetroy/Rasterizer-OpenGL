#include "triangle.h"

Coordinates::Coordinates(void)
{

    x = 0;
    y = 0;
    z = 0;
    s = 0;

    for(int i = 0; i < 3; i++){
        colors[i]  = 0;
        normals[i] = 0;
    }

}

Triangle::Triangle(void)
{
    for(int i = 0; i < 3; i++){

        X[i] = 0;
        Y[i] = 0;
        Z[i] = 0;
        S[i] = 0; 

        for(int k = 0; k < 3; k++){

            colors[i][k]  = 0;
            normals[i][k] = 0;
        }
    }

    triangletype.irregular  = false;
    triangletype.flattop    = false;
    triangletype.flatbottom = false;
    coordset                = false;

}

void
Triangle::Translate(Matrix totaltrans)
{

    double temp[4] = {0,0,0,0};

    for(int i = 0; i < 3; i++){
        double tcoordinates[4] = {X[i],Y[i],Z[i],1.};
        totaltrans.TransformPoint(tcoordinates, temp);
        X[i] = temp[0]/temp[3];
        Y[i] = temp[1]/temp[3];
        Z[i] = temp[2]/temp[3];
    }
}

void
Triangle::Print(ostream &o)
{
        char str[256];
        sprintf(str, "coord[0]:(%f,%f,%f) Shading:%f\n coord[1]:(%f,%f,%f) Shading:%f\n coord[2]:(%f,%f,%f) Shading:%f\n", coord[0].x,coord[0].y,coord[0].z,coord[0].s,coord[1].x,coord[1].y,coord[1].z,coord[1].s,coord[2].x,coord[2].y,coord[2].z,coord[2].s);
        o << str;
}


void
Triangle::InitializeTriangle(void)
{
    if((Y[0] == Y[1] && Y[1] == Y[2]) || (X[0] == X[1] && X[1] == X[2])){
        return;
    }

    if(Y[0] == Y[1] || Y[1] == Y[2] || Y[0] == Y[2]){
        Determinetype();
        SetCoordinates();
        return;
    }

    triangletype.irregular  = true;
    triangletype.flattop    = false;
    triangletype.flatbottom = false;
    SetCoordinates();

    return;
}

Triangle**
Triangle::SplitTriangle(void)
{
    if(!coordset){
        SetCoordinates();
    }

    Triangle** temparr  = new Triangle *[2];
    if(temparr == NULL)
        return NULL;

    Triangle *temp = new Triangle;
    if(temp == NULL)
        return NULL;

    Triangle *temp2 = new Triangle;
    if(temp == NULL)
        return NULL;

    double desiredy, desiredx, desiredz, desiredshading;
    double proportion;
    double desiredcolors[3];
    double desirednormals[3];
    desiredy = coord[1].y;

    proportion = (desiredy - coord[0].y) / (coord[2].y - coord[0].y);

    desiredx       = Lerp(coord[0].x, coord[2].x, proportion);
    desiredz       = Lerp(coord[0].z,coord[2].z,proportion);
    desiredshading = Lerp(coord[0].s,coord[2].s,proportion);

    for(int i = 0; i < 3; i ++){
        desiredcolors[i]   = Lerp(coord[0].colors[i], coord[2].colors[i], proportion);
        desirednormals[i]  = Lerp(coord[0].normals[i], coord[2].normals[i], proportion);
    }
    
    Normalize(desirednormals); 
    Normalize(coord[0].normals);
    Normalize(coord[1].normals);
    Normalize(coord[2].normals);

    temp->X[0] = coord[2].x;
    temp->Y[0] = coord[2].y;
    temp->Z[0] = coord[2].z;
    temp->S[0] = coord[2].s;
    memcpy(temp->colors[0],coord[2].colors,sizeof(double) * 3);
    memcpy(temp->normals[0],coord[2].normals,sizeof(double) * 3);
    temp->X[1] = coord[1].x;
    temp->Y[1] = coord[1].y;
    temp->Z[1] = coord[1].z;
    temp->S[1] = coord[1].s;
    memcpy(temp->colors[1],coord[1].colors,sizeof(double) * 3);
    memcpy(temp->normals[1],coord[1].normals,sizeof(double) * 3);
    temp->X[2] = desiredx;
    temp->Y[2] = desiredy;
    temp->Z[2] = desiredz;
    temp->S[2] = desiredshading;
    memcpy(temp->colors[2],desiredcolors,sizeof(double) * 3);
    memcpy(temp->normals[2],desirednormals,sizeof(double) * 3);

    temp2->X[0] = coord[1].x;
    temp2->Y[0] = coord[1].y;
    temp2->Z[0] = coord[1].z;
    temp2->S[0] = coord[1].s;
    memcpy(temp2->colors[0],coord[1].colors,sizeof(double) * 3);
    memcpy(temp2->normals[0],coord[1].normals,sizeof(double) * 3);
    temp2->X[1] = coord[0].x;
    temp2->Y[1] = coord[0].y;
    temp2->Z[1] = coord[0].z;
    temp2->S[1] = coord[0].s;
    memcpy(temp2->colors[1],coord[0].colors,sizeof(double) * 3);
    memcpy(temp2->normals[1],coord[0].normals,sizeof(double) * 3);
    temp2->X[2] = desiredx;
    temp2->Y[2] = desiredy;
    temp2->Z[2] = desiredz;
    temp2->S[2] = desiredshading;
    memcpy(temp2->colors[2],desiredcolors,sizeof(double) * 3);
    memcpy(temp2->normals[2],desirednormals,sizeof(double) * 3);

    temparr[0] = temp;
    temparr[1] = temp2;

    return temparr;
}

void
Triangle::Determinetype()
{
    if(Y[0] == Y[1]){
        if(Y[2] > Y[0]){
            triangletype.flatbottom = true;
        }else{
            triangletype.flattop = true;
        }
    }
    else if(Y[1] == Y[2]){
        if(Y[0] > Y[1]){
            triangletype.flatbottom = true;
        }else{
            triangletype.flattop = true;
        }
    }
    else if(Y[0] == Y[2]){
        if(Y[1] > Y[0]){
            triangletype.flatbottom = true;
        }else{
            triangletype.flattop = true;
        }
    }
    triangletype.irregular = false;
    return; 
}

void
Triangle::SetCoordinates(void)
{
    Coordinates left;
    Coordinates tip;
    Coordinates right;

    for(int i = 0; i < 3; i++){
        coord[i].x = X[i]; 
        coord[i].y = Y[i];
        coord[i].z = Z[i];
        coord[i].s = S[i];
        memcpy(coord[i].colors, colors[i], sizeof(double) * 3);
        memcpy(coord[i].normals, normals[i], sizeof(double) * 3);
    }

    if(triangletype.irregular){
        if((coord[0].y > coord[1].y && coord[0].y < coord[2].y) || (coord[0].y > coord[2].y && coord[0].y < coord[1].y)){ //if 0 is middle y
            tip.x = coord[0].x;
            tip.y = coord[0].y;
            tip.z = coord[0].z;
            tip.s = coord[0].s;
            memcpy(tip.colors,coord[0].colors, sizeof(double) * 3);
            memcpy(tip.normals,coord[0].normals, sizeof(double) * 3);
            if(coord[1].y < coord[2].y){
                left.x = coord[1].x;
                left.y = coord[1].y;
                left.z = coord[1].z;
                left.s = coord[1].s;
                memcpy(left.colors,coord[1].colors, sizeof(double) * 3);
                memcpy(left.normals,coord[1].normals, sizeof(double) * 3);
                right.x = coord[2].x;
                right.y = coord[2].y;
                right.z = coord[2].z;
                right.s = coord[2].s;
                memcpy(right.colors,coord[2].colors, sizeof(double) * 3);
                memcpy(right.normals,coord[2].normals, sizeof(double) * 3);
            } else {
                left.x = coord[2].x;
                left.y = coord[2].y;
                left.z = coord[2].z;
                left.s = coord[2].s;
                memcpy(left.colors,coord[2].colors, sizeof(double) * 3);
                memcpy(left.normals,coord[2].normals, sizeof(double) * 3);
                right.x = coord[1].x;
                right.y = coord[1].y;
                right.z = coord[1].z;
                right.s = coord[1].s;
                memcpy(right.colors,coord[1].colors, sizeof(double) * 3);
                memcpy(right.normals,coord[1].normals, sizeof(double) * 3);
            }
        }
        else if((coord[1].y > coord[0].y && coord[1].y < coord[2].y) || (coord[1].y > coord[2].y && coord[1].y < coord[0].y)){ //if 1 is middle y
            tip.x = coord[1].x;
            tip.y = coord[1].y;
            tip.z = coord[1].z;
            tip.s = coord[1].s;
            memcpy(tip.colors,coord[1].colors, sizeof(double) * 3);
            memcpy(tip.normals,coord[1].normals, sizeof(double) * 3);
            if(coord[0].y < coord[2].y){
                left.x = coord[0].x;
                left.y = coord[0].y;
                left.z = coord[0].z;
                left.s = coord[0].s;
                memcpy(left.colors,coord[0].colors, sizeof(double) * 3);
                memcpy(left.normals,coord[0].normals, sizeof(double) * 3);
                right.x = coord[2].x;
                right.y = coord[2].y;
                right.z = coord[2].z;
                right.s = coord[2].s;
                memcpy(right.colors,coord[2].colors, sizeof(double) * 3);
                memcpy(right.normals,coord[2].normals, sizeof(double) * 3);
            } else {
                left.x = coord[2].x;
                left.y = coord[2].y;
                left.z = coord[2].z;
                left.s = coord[2].s;
                memcpy(left.colors,coord[2].colors, sizeof(double) * 3);
                memcpy(left.normals,coord[2].normals, sizeof(double) * 3);
                right.x = coord[0].x;
                right.y = coord[0].y;
                right.z = coord[0].z;
                right.s = coord[0].s;
                memcpy(right.colors,coord[0].colors, sizeof(double) * 3);
                memcpy(right.normals,coord[0].normals, sizeof(double) * 3);
            }
        }       
        else if((coord[2].y > coord[1].y && coord[2].y < coord[0].y) || (coord[2].y > coord[0].y && coord[2].y < coord[1].y)){ //if 2 is middle y
            tip.x = coord[2].x;
            tip.y = coord[2].y;
            tip.z = coord[2].z;
            tip.s = coord[2].s;
            memcpy(tip.colors,coord[2].colors, sizeof(double) * 3);
            memcpy(tip.normals,coord[2].normals, sizeof(double) * 3);
            if(coord[1].y < coord[0].y){
                left.x  = coord[1].x;
                left.y  = coord[1].y;
                left.z  = coord[1].z;
                left.s = coord[1].s;
                memcpy(left.colors,coord[1].colors, sizeof(double) * 3);
                memcpy(left.normals,coord[1].normals, sizeof(double) * 3);
                right.x = coord[0].x;
                right.y = coord[0].y;
                right.z = coord[0].z;
                right.s = coord[0].s;
                memcpy(right.colors,coord[0].colors, sizeof(double) * 3);
                memcpy(right.normals,coord[0].normals, sizeof(double) * 3);
            } else {
                left.x  = coord[0].x;
                left.y  = coord[0].y;
                left.z  = coord[0].z;
                left.s = coord[0].s;
                memcpy(left.colors,coord[0].colors, sizeof(double) * 3);
                memcpy(left.normals,coord[0].normals, sizeof(double) * 3);
                right.x = coord[1].x;
                right.y = coord[1].y;
                right.z = coord[1].z;
                right.s = coord[1].s;
                memcpy(right.colors,coord[1].colors, sizeof(double) * 3);
                memcpy(right.normals,coord[1].normals, sizeof(double) * 3);
            }
        }
    }

    if(triangletype.flattop || triangletype.flatbottom){
        if(coord[0].y == coord[1].y){
            tip.x = coord[2].x;
            tip.y = coord[2].y;
            tip.z = coord[2].z;
            tip.s = coord[2].s;
            memcpy(tip.colors,coord[2].colors, sizeof(double) * 3);
            memcpy(tip.normals,coord[2].normals, sizeof(double) * 3);
            if(coord[0].x > coord[1].x){
                left.x = coord[1].x;
                left.y = coord[1].y;
                left.z = coord[1].z;
                left.s = coord[1].s;
                memcpy(left.colors,coord[1].colors, sizeof(double) * 3);
                memcpy(left.normals,coord[1].normals, sizeof(double) * 3);
                right.x = coord[0].x;
                right.y = coord[0].y;
                right.z = coord[0].z;
                right.s = coord[0].s;
                memcpy(right.colors,coord[0].colors, sizeof(double) * 3);
                memcpy(right.normals,coord[0].normals, sizeof(double) * 3);
            } else {
                left.x = coord[0].x;
                left.y = coord[0].y;
                left.z = coord[0].z;
                left.s = coord[0].s;
                memcpy(left.colors,coord[0].colors, sizeof(double) * 3);
                memcpy(left.normals,coord[0].normals, sizeof(double) * 3);
                right.x = coord[1].x;
                right.y = coord[1].y;
                right.z = coord[1].z;
                right.s = coord[1].s;
                memcpy(right.colors,coord[1].colors, sizeof(double) * 3);
                memcpy(right.normals,coord[1].normals, sizeof(double) * 3);
            }
        }
        else if(coord[1].y == coord[2].y){
            tip.x = coord[0].x;
            tip.y = coord[0].y;
            tip.z = coord[0].z;
            tip.s = coord[0].s;
            memcpy(tip.colors,coord[0].colors, sizeof(double) * 3);
            memcpy(tip.normals,coord[0].normals, sizeof(double) * 3);
            if(coord[1].x > coord[2].x){
                left.x = coord[2].x;
                left.y = coord[2].y;
                left.z = coord[2].z;
                left.s = coord[2].s;
                memcpy(left.colors,coord[2].colors, sizeof(double) * 3);
                memcpy(left.normals,coord[2].normals, sizeof(double) * 3);
                right.x = coord[1].x;
                right.y = coord[1].y;
                right.z = coord[1].z;
                right.s = coord[1].s;
                memcpy(right.colors,coord[1].colors, sizeof(double) * 3);
                memcpy(right.normals,coord[1].normals, sizeof(double) * 3);
            } else {
                left.x = coord[1].x;
                left.y = coord[1].y;
                left.z = coord[1].z;
                left.s = coord[1].s;
                memcpy(left.colors,coord[1].colors, sizeof(double) * 3);
                memcpy(left.normals,coord[1].normals, sizeof(double) * 3);
                right.x = coord[2].x;
                right.y = coord[2].y;
                right.z = coord[2].z;
                right.s = coord[2].s;
                memcpy(right.colors,coord[2].colors, sizeof(double) * 3);
                memcpy(right.normals,coord[2].normals, sizeof(double) * 3);
            }
        }
        else if(coord[0].y == coord[2].y){
            tip.x = coord[1].x;
            tip.y = coord[1].y;
            tip.z = coord[1].z;
            tip.s = coord[1].s;
            memcpy(tip.colors,coord[1].colors, sizeof(double) * 3);
            memcpy(tip.normals,coord[1].normals, sizeof(double) * 3);
            if(coord[0].x > coord[2].x){
                left.x = coord[2].x;
                left.y = coord[2].y;
                left.z = coord[2].z;
                left.s = coord[2].s;
                memcpy(left.colors,coord[2].colors, sizeof(double) * 3);
                memcpy(left.normals,coord[2].normals, sizeof(double) * 3);
                right.x = coord[0].x;
                right.y = coord[0].y;
                right.z = coord[0].z;
                right.s = coord[0].s;
                memcpy(right.colors,coord[0].colors, sizeof(double) * 3);
                memcpy(right.normals,coord[0].normals, sizeof(double) * 3);
            } else {
                left.x = coord[0].x;
                left.y = coord[0].y;
                left.z = coord[0].z;
                left.s = coord[0].s;
                memcpy(left.colors,coord[0].colors, sizeof(double) * 3);
                memcpy(left.normals,coord[0].normals, sizeof(double) * 3);
                right.x = coord[2].x;
                right.y = coord[2].y;
                right.z = coord[2].z;
                right.s = coord[2].s;
                memcpy(right.colors,coord[2].colors, sizeof(double) * 3);
                memcpy(right.normals,coord[2].normals, sizeof(double) * 3);
            }
        }
    }

    coord[0] = left;
    coord[1] = tip;
    coord[2] = right;
    coordset = true;

    return;
}