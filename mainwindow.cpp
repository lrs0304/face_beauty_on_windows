#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QDebug>
#include <QFileDialog>
#include <QSlider>
#include <chrono>
#include <omp.h>
#include <math.h>

void MeanCovMapCalculate(float*data, int width, int height, float *corrIP, int radius);
void f_EPMFilter(const uchar* srcData, uchar* dstData, int nWidth, int nHeight, int nStride, int radius, float delta, float light);

//参考 https://blog.csdn.net/Trent1985/article/details/80802144#commentsedit
MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    initWidgets();
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::initWidgets()
{
    connect(ui->pushButton, &QPushButton::clicked, this, [=] {

        //QFileDialog* loadFileDialog = new QFileDialog(this);
        //loadFileDialog->setWindowTitle(QString("Select An Picture"));
        //loadFileDialog->setDirectory("E:/Rison/Pictures");
        ////设置可以选择多个文件,默认为只能选择一个文件QFileDialog::ExistingFiles
        ////fileDialog->setFileMode(QFileDialog::ExistingFiles);
        //loadFileDialog->setNameFilter(tr("Images(*.png *.jpg *.jpeg *.bmp)"));
        //if(loadFileDialog->exec())
        {
            //auto fileNames = loadFileDialog->selectedFiles();
            srcImg = QImage("E:/Rison/Pictures/fff.png");
            outImg = srcImg.copy();
            width = srcImg.width();
            height = srcImg.height();
            qDebug() << "byte count " << srcImg.byteCount();
            ui->input->setPixmap(QPixmap::fromImage(srcImg));
        }
        //loadFileDialog->deleteLater();
    });

    QSlider *slider = new QSlider(this);
    slider->setGeometry(320, 270, 200, 20);
    slider->setOrientation(Qt::Horizontal);
    slider->setMinimum(15);
    slider->setMaximum(50);
    connect(slider, &QSlider::valueChanged, this, [=](int value) {

        ///////////////////////////////////////////////
        auto beginClock = std::chrono::system_clock::now();
        f_EPMFilter(srcImg.bits(), outImg.bits(), width, height, width, 10, 0.001f,value/10.0f);
        auto endClock = std::chrono::system_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(endClock - beginClock);
        qDebug("cost time: %.3lf mms", duration.count() / 1000.0f);//输出图像处理花费时间信息
        ///////////////////////////////////////////////

        ui->output->setPixmap(QPixmap::fromImage(outImg));
    });
}

//使用积分图用于快速计算均值，进而计算方差
//参考 https://blog.csdn.net/fb_help/article/details/88534427
void MeanCovMapCalculate(float*data, int width, int height, float *corrIP, int radius)
{
    //sump需要比原来的多1，因为用(1,1)处的值表示(0,0)的值，(x+1,y+1)表示(0,0)~(x,y)直接的总和
    //qDebug() <<"--->"<<__FUNCTION__ <<":"<<__LINE__ <<" size: "<<width<<"x"<<height;
    int mapw=width+1;
    int maph=height+1;
    float *summap = (float*)malloc(sizeof(float) * mapw * maph);
    memset(summap, 0, sizeof(float) * mapw * maph);
    //qDebug() <<"--->"<<__FUNCTION__ <<":"<<__LINE__;

    for (int x=0; x<width; x++)
    {
        for (int y=0; y<height; y++)
        {
            int index = mapw*y+x;
            //           x+1,y+1                  x+1,y                  x,y+1               x,y              x,y
            //summap[mapw*(x+1)+y+1] = summap[mapw*(x+1)+y] + summap[mapw*x+y+1] - summap[mapw*x+y] + data[width*x+y];
            //qDebug() <<  "--->x"<<(x+1)<<",y"<<(y+1) <<"|"<< summap[mapw*(x+1)+y+1]  <<  summap[mapw*(x+1)+y]  <<  summap[mapw*x+y+1]  <<  summap[mapw*x+y]  <<  data[width*x+y];
            summap[index+mapw+1] = summap[index+mapw] + summap[index+1] - summap[index] + data[index-y];
        }
    }

    //qDebug() <<"--->"<<__FUNCTION__ <<":"<<__LINE__;
    int left=radius/2;
    int right=radius-left;
    int top = radius/2;
    int bottom = radius-top;
    for (int x=0; x<width; ++x) {
        for (int y=0; y<height; ++y) {
            int x1 = x-left;
            if(x1<0){x1=0;}

            int x2 = x+right;
            if(x2>width){x2=width;}

            int y1 = y-top;
            if(y1<0){y1=0;}

            int y2=y+bottom;
            if(y2>height){y2=height;}

            int n=(x2-x1+1)*(y2-y1+1);
            //qDebug() <<"--->" <<(mapw*y2+x2)<<",total:"<<(mapw * maph);
            corrIP[width*y+x] = (summap[mapw*y2+x2]-summap[mapw*y2+x1]-summap[mapw*y1+x2]+summap[mapw*y1+x1])/n;
        }
    }

    //qDebug() <<"--->"<<__FUNCTION__ <<":"<<__LINE__;
    free(summap);
    //qDebug() <<"--->"<<__FUNCTION__ <<":"<<__LINE__;
}

unsigned char inline CLIP3(int v, int min, int max)
{
    if (v < min) { return min; }
    else if (v > max) { return max; }
    else return v;
}

float inline light1(float valueDividedBy255, float beta)
{
    return log(valueDividedBy255*(beta - 1) + 1) / log(beta);
}

//Edge Preserved mean filter
int EPMFilter(unsigned char* srcData, int width, int height, int radius, float delta, float light)
{
    int size = width*height;
    float *data = (float*)malloc(sizeof(float) * size);
    float *meanI = (float*)malloc(sizeof(float) * size);
    float *corrI = (float*)malloc(sizeof(float) * size);

    //使用积分图计算平均值
    for (int i = 0; i < size; i++)
    {
        data[i] = (float)srcData[i] / 255.0f;
    }
    MeanCovMapCalculate(data, width, height, meanI, radius);

    //曲线美白
    for (int i = 0; i < size; i++)
    {
        srcData[i] = (uchar)CLIP3(light1(data[i], light) * 255.0f, 0, 255);
    }

    //计算方差
    for (int i = 0; i < size; i++)
    {
        data[i] *= data[i];
    }
    MeanCovMapCalculate(data, width, height, corrI, radius);
    for (int i = 0; i < size; i++)
    {
        corrI[i] = corrI[i] - meanI[i] * meanI[i];
    }

    //均方差磨皮，暂存到highPass里面，等待处理
    for (int i = 0; i < size; i++)
    {
        float t = meanI[i] + (corrI[i] * (srcData[i] / 255.0f - meanI[i]) / (corrI[i] + delta));
        srcData[i] = (uchar)CLIP3(t * 255.0f, 0, 255);
    }

    free(data);
    free(meanI);
    free(corrI);
    return 0;
}

//3通道处理
void f_EPMFilter(const uchar* srcData, uchar* dstData, int nWidth, int nHeight, int nStride, int radius, float delta, float light)
{
    Q_UNUSED(nStride);
    if (srcData == nullptr) { return; }

    int64_t rBytesCount = sizeof(uchar) * nWidth * nHeight;
    uchar* rData = (uchar*)malloc(rBytesCount);
    uchar* gData = (uchar*)malloc(rBytesCount);
    uchar* bData = (uchar*)malloc(rBytesCount);
    const uchar* pSrc = srcData;

    uchar* pR = rData;
    uchar* pG = gData;
    uchar* pB = bData;
    for (int y = 0; y < nHeight; y++)
    {
        for (int x = 0; x < nWidth; x++)
        {
            *pR = pSrc[0];
            *pG = pSrc[1];
            *pB = pSrc[2];
            pR++;
            pG++;
            pB++;
            pSrc += 4;
        }
    }

    //并行处理
#pragma omp parallel sections
    {
#pragma omp section
        EPMFilter(rData, nWidth, nHeight, radius, delta, light);
#pragma omp section
        EPMFilter(gData, nWidth, nHeight, radius, delta, light);
#pragma omp section
        EPMFilter(bData, nWidth, nHeight, radius, delta, light);
    }

    uchar* pOut = dstData;
    pR = rData;
    pG = gData;
    pB = bData;

    for (int j = 0; j < nHeight; j++)
    {
        for (int i = 0; i < nWidth; i++)
        {
            pOut[0] = *pR;
            pOut[1] = *pG;
            pOut[2] = *pB;
            pR++;
            pG++;
            pB++;
            pOut += 4;
        }
    }

    free(rData);
    free(gData);
    free(bData);
}

/*
void inline RGBToYCbCr(unsigned char R, unsigned char G, unsigned char B, int &Y, int &Cb, int &Cr)
{
    Y = (unsigned char)(0.257*R + 0.564*G + 0.098*B + 16);
    Cb = (unsigned char)(-0.148*R - 0.291*G + 0.439*B + 128);
    Cr = (unsigned char)(0.439*R - 0.368*G + 0.071*B + 128);
}

void inline YCbCrToRGB(unsigned char Y, unsigned char Cb, unsigned char Cr, int &R, int &G, int &B)
{
    R = 1.164*(Y - 16) + 1.596*(Cr - 128);
    G = 1.164*(Y - 16) - 0.392*(Cb - 128) - 0.813*(Cr - 128);
    B = 1.164*(Y - 16) + 2.017*(Cb - 128);
}

//只处理Y通道
void f_EPMFilterOneChannel(unsigned char* srcData, int nWidth, int nHeight, int nStride, int radius, float delta)
{
    if (srcData == NULL)
    {
        return;
    }
    unsigned char* yData = (unsigned char*)malloc(sizeof(unsigned char) * nWidth * nHeight);
    unsigned char* cbData = (unsigned char*)malloc(sizeof(unsigned char) * nWidth * nHeight);
    unsigned char* crData = (unsigned char*)malloc(sizeof(unsigned char) * nWidth * nHeight);
    unsigned char* pSrc = srcData;
    int Y, CB, CR;
    unsigned char* pY = yData;
    unsigned char* pCb = cbData;
    unsigned char* pCr = crData;
    for(int j = 0; j < nHeight; j++)
    {
        for(int i = 0; i < nWidth; i++)
        {
            RGBToYCbCr(pSrc[2],pSrc[1],pSrc[0], Y, CB, CR);
            *pY = Y;
            *pCb = CB;
            *pCr = CR;
            pY++;
            pCb++;
            pCr++;
            pSrc += 4;
        }
    }
    EPMFilter(yData, nWidth, nHeight, radius, delta);
    pSrc = srcData;
    pY = yData;
    pCb = cbData;
    pCr = crData;
    int R, G, B;
    for(int j = 0; j < nHeight; j++)
    {
        for(int i = 0; i < nWidth; i++)
        {
            YCbCrToRGB(*pY, *pCb, *pCr, R, G, B);
            pSrc[0] = B;
            pSrc[1] = G;
            pSrc[2] = R;
            pY++;
            pCb++;
            pCr++;
            pSrc += 4;
        }
    }free(yData);free(cbData);free(crData);
}
*/
