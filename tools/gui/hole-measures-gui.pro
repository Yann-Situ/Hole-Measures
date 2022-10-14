#-------------------------------------------------
#
# Project created by QtCreator 2018-10-02T16:25:45
#
#-------------------------------------------------

QT       += core gui
QT       += opengl
QT       += charts

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = hole-measures-gui
TEMPLATE = app

# The following define makes your compiler emit warnings if you use
# any feature of Qt which has been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0


unix:!macx {
    LIBS += -lglut -lGLU
    LIBS += -L$$PWD/OpenMesh/liblinux/ -lOpenMeshCore

    INCLUDEPATH += $$PWD/OpenMesh/inc/
    DEPENDPATH += $$PWD/OpenMesh/inc/
    DEPENDPATH += $$PWD/OpenMesh/liblinux/
}

macx: {
    INCLUDEPATH += $$PWD/../OpenMesh/inc/
    LIBS += -L$$PWD/../OpenMesh/libosx/ -lOpenMeshCore -lOpenMeshTools
}

SOURCES += \
    chartview.cpp \
        main.cpp \
        mainwindow.cpp \
        meshviewerwidget.cpp \

HEADERS += \
    chartview.h \
        mainwindow.h \
        meshviewerwidget.h \

FORMS += \
        mainwindow.ui


