#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QImage>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private:
    void initWidgets();
    Ui::MainWindow *ui;
    QImage srcImg;
    QImage outImg;
    int width, height;
};

#endif // MAINWINDOW_H
