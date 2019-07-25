#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QDebug>
#include <QFileDialog>

void MeanCovMapCalculate(float*data, int width, int height, float *corrIP, int radius);
void f_EPMFilter(uchar* srcData, int nWidth, int nHeight, int nStride, int radius, float delta);

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    initWidgets();
    float *hehe = (float*)malloc(sizeof(float) * 5 * 5);
    for(int i=0;i<5;i++)
    {
        for(int j=0;j<5;j++)
        {
            hehe[i*5+j]=1;
        }
    }
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::initWidgets()
{
    connect(ui->pushButton,&QPushButton::clicked, this, [=]{

        QFileDialog* loadFileDialog = new QFileDialog(this);
        loadFileDialog->setWindowTitle(QString("请选择图片"));
        loadFileDialog->setDirectory("E:/Rison/Desktop");
        //设置可以选择多个文件,默认为只能选择一个文件QFileDialog::ExistingFiles
        //fileDialog->setFileMode(QFileDialog::ExistingFiles);
        loadFileDialog->setNameFilter(tr("Images(*.png *.jpg *.jpeg *.bmp)"));
        if(loadFileDialog->exec())
        {
            auto fileNames = loadFileDialog->selectedFiles();
            QPixmap pixmap(fileNames.first());
            ui->input->setPixmap(pixmap);//(pixmap.scaled(ui->input->size(), Qt::KeepAspectRatio));
        }

        loadFileDialog->deleteLater();
    });

    connect(ui->pushButton_2, &QPushButton::clicked, this, [=]{

        QImage img=QImage("E:/Rison/Desktop/ddd.png").convertToFormat(QImage::Format_ARGB32);// = ui->input->pixmap()->toImage().convertToFormat(QImage::Format_ARGB32);


        ///////////////////////////////////////////////
        qDebug() <<"byte count " << img.byteCount();
        //uchar* data = (uchar*)malloc(sizeof(uchar) * img.byteCount());;
        //for(int y=0;y<img.height();y++)
        //{
        //    uchar*tmp = img.scanLine(y);
        //    memcpy(data+img.bytesPerLine()*y, tmp, img.bytesPerLine());
        //}
        //
        //QImage hehe(img.bits(),img.width(),img.height(),QImage::Format_ARGB32);

        f_EPMFilter(img.bits(), img.width(), img.height(), img.width(), 10, 0.04f);
        ui->output->setPixmap(QPixmap::fromImage(img));
        //free(data);
        //img.save("d:/rison.jpg");

    });
}

//使用积分图用于快速计算均值，进而计算方差
void MeanCovMapCalculate(float*data, int width, int height, float *corrIP, int radius)
{
    //sump需要比原来的多1，因为用(1,1)处的值表示(0,0)的值，(x+1,y+1)表示(0,0)~(x,y)直接的总和
    qDebug() <<"--->"<<__FUNCTION__ <<":"<<__LINE__ <<" size: "<<width<<"x"<<height;
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
    qDebug() <<"--->"<<__FUNCTION__ <<":"<<__LINE__;
}

unsigned char inline CLIP3(int v, int min, int max)
{
    if(v<min){return v;}
    else if(v>max){return max;}
    else return v;
}

//Edge Preserved mean filter
int EPMFilter(unsigned char* srcData, int width ,int height, int radius, float delta)
{
    float *data = (float*)malloc(sizeof(float) * width * height);
    float *meanIP = (float*)malloc(sizeof(float) * width * height);
    float *corrIP = (float*)malloc(sizeof(float) * width * height);
    for(int i = 0; i < width * height; i++)
    {
        data[i] = (float)srcData[i] / 255.0f;
    }
    //mean and cov compute
    MeanCovMapCalculate(data, width, height, meanIP, radius);
    for(int i = 0; i < width * height; i++)
    {
        data[i] *= data[i];
    }
    //mean and cov compute
    MeanCovMapCalculate(data, width, height, corrIP, radius);
    for(int i = 0; i < width * height; i++)
    {
        corrIP[i] = corrIP[i] - meanIP[i] * meanIP[i];
    }
    for(int i = 0; i < width * height; i++)
    {
        float t = meanIP[i] + (corrIP[i] * (srcData[i] / 255.0f - meanIP[i]) / (corrIP[i] + delta));
        srcData[i] = (unsigned char)(CLIP3(t * 255.0f, 0, 255));
    }
    free(data);
    free(meanIP);
    free(corrIP);
    return 0;
}

//4通道处理
void f_EPMFilter(uchar* srcData, int nWidth, int nHeight, int nStride, int radius, float delta)
{
    Q_UNUSED(nStride);

    if (srcData == NULL) { return; }

    qDebug() <<"--->"<<__FUNCTION__ <<":"<<__LINE__;
    uchar* rData = (uchar*)malloc(sizeof(uchar) * nWidth * nHeight);
    uchar* gData = (uchar*)malloc(sizeof(uchar) * nWidth * nHeight);
    uchar* bData = (uchar*)malloc(sizeof(uchar) * nWidth * nHeight);
    uchar* pSrc = srcData;

    qDebug() <<"---->width:"<<nWidth<<",height:"<<nHeight;
    uchar* pR = rData;
    uchar* pG = gData;
    uchar* pB = bData;
    qDebug() <<"--->"<<__FUNCTION__ <<":"<<__LINE__;
    for(int y = 0; y < nHeight; y++)
    {
        for(int x = 0; x < nWidth; x++)
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

    qDebug() <<"--->"<<__FUNCTION__ <<":"<<__LINE__;
    //并行处理
    //#pragma omp parallel sections  num_threads(omp_get_num_procs())
    {
    //    #pragma omp  section
            EPMFilter(rData, nWidth, nHeight, radius, delta);
    //    #pragma omp  section
            EPMFilter(gData, nWidth, nHeight, radius, delta);
    //    #pragma omp  section
            EPMFilter(bData, nWidth, nHeight, radius, delta);
    }
    qDebug() <<"--->"<<__FUNCTION__ <<":"<<__LINE__;
    for(int i=0;i<15;i++)
    {
        qDebug() <<"~~~~~~"<<srcData[i*3] <<" vs " << rData[i];
    }
    pSrc = srcData;
    pR = rData;
    pG = gData;
    pB = bData;
    qDebug() <<"--->"<<__FUNCTION__ <<":"<<__LINE__;

    for(int j = 0; j < nHeight; j++)
    {
        for(int i = 0; i < nWidth; i++)
        {
            pSrc[0] = *pR;
            pSrc[1] = *pG;
            pSrc[2] = *pB;
            pR++;
            pG++;
            pB++;
            pSrc += 4;
        }
    }
    qDebug() <<"--->"<<__FUNCTION__ <<":"<<__LINE__;
    free(rData);
    free(gData);
    free(bData);
    qDebug() <<"--->"<<__FUNCTION__ <<":"<<__LINE__;
}

void inline RGBToYCbCr(unsigned char R,unsigned char G,unsigned char B, int &Y, int &Cb, int &Cr)
{
    Y =  (unsigned char)( 0.257*R+0.564*G+0.098*B+16);
    Cb = (unsigned char)(-0.148*R-0.291*G+0.439*B+128);
    Cr = (unsigned char)( 0.439*R-0.368*G+0.071*B+128);
}

void inline YCbCrToRGB(unsigned char Y,unsigned char Cb,unsigned char Cr, int &R, int &G, int &B)
{
    R = 1.164*(Y-16)+1.596*(Cr-128);
    G = 1.164*(Y-16)-0.392*(Cb-128)-0.813*(Cr-128);
    B = 1.164*(Y-16)+2.017*(Cb-128);
}
/*
//单通道处理
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
