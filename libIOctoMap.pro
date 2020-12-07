#-------------------------------------------------
#
# Project created by QtCreator 2018-07-03T08:38:32
#
#-------------------------------------------------
debug_and_release {
    CONFIG -= debug_and_release
    CONFIG += debug_and_release
}

# ensure one "debug" or "release" in CONFIG so they can be used as
# conditionals instead of writing "CONFIG(debug, debug|release)"...
CONFIG(debug, debug|release) {
    CONFIG -= debug release
    CONFIG += debug
    }
CONFIG(release, debug|release) {
        CONFIG -= debug release
        CONFIG += release
}

#QT       += opengl
QT += widgets

TARGET = libIOctoMap
TEMPLATE = lib

DEFINES += LIBIOCTOMAP_LIBRARY

SOURCES += libIOctoMap.cpp

HEADERS += libIOctoMap.h\
        libioctomap_global.h

OCTOMAP_PATH=E:/Librerias/octomap-1.9.0/install

INCLUDEPATH += $$OCTOMAP_PATH/include

DESTDIR_RELEASE=C:\dev\release
DESTDIR_DEBUG=C:\dev\debug
#DESTDIR_RELEASE=C:\dev\release\PCL
#DESTDIR_DEBUG=C:\dev\debug\PCL
#OSGEO4W_PATH="C:\Program Files\QGIS Essen"
debug{
    DESTDIR = $$DESTDIR_DEBUG
    LIBS += -L$$DESTDIR_DEBUG
    LIBS += $$OCTOMAP_PATH\debug\lib\octomap.lib
    LIBS += $$OCTOMAP_PATH\debug\lib\octomath.lib
}else{
    DESTDIR = $$DESTDIR_RELEASE
    LIBS += -L$$DESTDIR_RELEASE
    LIBS += $$OCTOMAP_PATH\release\lib\octomap.lib
    LIBS += $$OCTOMAP_PATH\release\lib\octomath.lib
}

#unix {
#    target.path = /usr/lib
#    INSTALLS += target
#}
