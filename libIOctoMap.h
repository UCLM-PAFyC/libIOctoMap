#ifndef LIBIOCTOMAP_H
#define LIBIOCTOMAP_H

#include "libioctomap_global.h"

#include <QMap>
#include <QVector>
#include <QString>
#include <QWidget>

#define LIBIOCTOMAP_OPENVIS                     "octovis/octovis.exe"
#define LIBIOCTOMAP_BT_FILE_EXTENSION           ".bt"
#define LIBIOCTOMAP_OT_FILE_EXTENSION           ".ot"
#define LIBIOCTOMAP_NUMBER_OF_POINTS_TO_PROCESS_BY_STEP     10000

#define LIBIOCTOMAP_LINEAR_TOLERANCE                        0.0001

#include <math.h>
#include <stdlib.h>

// this is mimicing gtest expressions

#define EXPECT_TRUE(args) {                                             \
    if (!(args)) { fprintf(stderr, "test failed (EXPECT_TRUE) in %s, line %d\n", __FILE__, __LINE__); \
      exit(1);                                                         \
    } }

#define EXPECT_FALSE(args) {                                             \
    if (args) { fprintf(stderr, "test failed (EXPECT_FALSE) in %s, line %d\n", __FILE__, __LINE__); \
      exit(1);                                                         \
    } }

#define EXPECT_EQ(a,b) {                                                \
    if (!(a == b)) { std::cerr << "test failed: " <<a<<"!="<<b<< " in " \
                      << __FILE__ << ", line " <<__LINE__ << std::endl; \
      exit(1);                                                          \
    } }

#define EXPECT_FLOAT_EQ(a,b) {                                          \
    if (!(fabs(a-b) <= 1e-5)) { fprintf(stderr, "test failed: %f != %f in %s, line %d\n", a, b, __FILE__, __LINE__); \
      exit(1);                                                         \
    } }

#define EXPECT_NEAR(a,b,prec) {                                         \
    if (!(fabs(a-b) <= prec)) { fprintf(stderr, "test failed: |%f - %f| > %f in %s, line %d\n", a, b, prec, __FILE__, __LINE__); \
      exit(1);                                                         \
    } }

class LIBIOCTOMAPSHARED_EXPORT libIOctoMap
{

public:
    inline libIOctoMap() {;};
    static inline libIOctoMap * getInstance(void )
    {
        if ( mInstance == 0 ) mInstance = new libIOctoMap;
            return mInstance;
    };
    bool createOctoMap(QVector<QVector<double> >& coordinates,
                       QVector<QVector<unsigned short> >& colors,
                       QVector<unsigned short>& sourceIds,
                       QMap<int,QVector<double> >& scanPositions,
                       float resolution,
                       bool computeFreeVoxels,
                       double& volume,
                       bool computeMetrics,
                       bool openViewer,
                       QString outputFileName,
                       QString& strError,
                       QWidget *ptrWidget=NULL);
    bool openFileInViewer(QString fileName,
                    QString& strError);
private:
    double max(double v1,double v2);
    double min(double v1,double v2);
    bool getNewOriginForRayInBoundingBox(double xs,double ys,double zs, // start
                                         double xe,double ye,double ze, // end
                                         double minX, double minY,
                                         double maxX, double maxY,
                                         double& xo,
                                         double& yo,
                                         double& zo);
    static libIOctoMap * mInstance;
};

#endif // LIBIOCTOMAP_H
