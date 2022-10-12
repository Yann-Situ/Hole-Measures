#include "meshviewerwidget.h"

MeshViewerWidget::MeshViewerWidget(QWidget*_parent) : QGLWidget(_parent)
{
    triToDraw = 0;
    linesToDraw = 0;
    pointsToDraw = 0;

    setMouseTracking(true);
    setFocus();
}

MeshViewerWidget::MeshViewerWidget( QGLFormat& _fmt, QWidget* _parent ) : QGLWidget( _fmt, _parent )
{
    setMouseTracking(true);
    setFocus();
}

void MeshViewerWidget::initializeGL()
{
    glClearColor(1.0, 1.0, 1.0, 1.0);
    glEnable( GL_DEPTH_TEST );

    glLoadIdentity();

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glGetDoublev(GL_MODELVIEW_MATRIX, modelview_matrix_);
    set_scene_pos(Vec3f(0.0, 0.0, 0.0), 1.0);

    glEnable( GL_MULTISAMPLE );

    // glShadeModel(GL_SMOOTH);//Enable shadow smoothing
    //
    // GLfloat LightAmbient[] = {0.5f, 0.5f, 0.6f, 1.0f};//Ambient light parameters
    // GLfloat LightDiffuse[] = {0.6f, 0.6f, 0.5f, 1.0f};//Diffuse light parameters
    // GLfloat LightPosition[] = {0.0f, 0.0f, 1.0f, 0.0f};//light source position
    // glLightfv(GL_LIGHT1, GL_AMBIENT, LightAmbient);//Set the ambient light
    // glLightfv(GL_LIGHT1, GL_DIFFUSE, LightDiffuse);//Set diffuse light
    // glLightfv(GL_LIGHT1, GL_POSITION, LightPosition);//Set the light source position
    // glEnable(GL_LIGHT1);//Start light source No. 1
    //
    // glEnable(GL_LIGHTING);
    // glEnable(GL_COLOR_MATERIAL);
    // glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
}

void MeshViewerWidget::translate( const OpenMesh::Vec3f& _trans )
{
    glLoadIdentity();
    glTranslated( _trans[0], _trans[1], _trans[2] );
    glMultMatrixd( modelview_matrix_ );
    glGetDoublev( GL_MODELVIEW_MATRIX, modelview_matrix_);
}

void MeshViewerWidget::resizeGL( int _w, int _h )
{
    update_projection_matrix();
    glViewport(0, 0, _w, _h);
    updateGL();
}

void MeshViewerWidget::reloadPOV()
{
    glLoadIdentity();
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glGetDoublev(GL_MODELVIEW_MATRIX, modelview_matrix_);
    set_scene_pos(Vec3f(0.0, 0.0, 0.0), 1.0);
    updateGL();
}

void MeshViewerWidget::loadMesh(GLfloat* verts, GLfloat* colors, int nVerts, GLuint* triangles, int nTriangles)
{
    GLfloat* vertsColsArray = new GLfloat[nVerts * 2];

    for(int i = 0; i < nVerts; i = i + 3)
    {
        int j = 2 * i;
        vertsColsArray[j] = colors[i] / 255.0;
        vertsColsArray[j+1] = colors[i+1] / 255.0;
        vertsColsArray[j+2] = colors[i+2] / 255.0;

        vertsColsArray[j+3] = verts[i];
        vertsColsArray[j+4] = verts[i+1];
        vertsColsArray[j+5] = verts[i+2];
    }

    glGenBuffers( 2, TriDataBuffers );

    glBindBuffer(GL_ARRAY_BUFFER, TriDataBuffers[0]);
    glBufferData(GL_ARRAY_BUFFER, nVerts * 2 * sizeof(GLfloat), vertsColsArray, GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, TriDataBuffers[1]);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, nTriangles * sizeof(GLuint), triangles, GL_STATIC_DRAW);

    triToDraw = nTriangles;

    init = true;

    delete[] vertsColsArray;

    updateGL();
}

void MeshViewerWidget::loadBalls(GLfloat *verts, GLfloat *colors, int nVerts, GLuint *triangles, int nTriangles)
{
    GLfloat* vertsColsArray = new GLfloat[nVerts * 2];

    for(int i = 0; i < nVerts; i = i + 3)
    {
        int j = 2 * i;
        vertsColsArray[j] = colors[i] / 255.0;
        vertsColsArray[j+1] = colors[i+1] / 255.0;
        vertsColsArray[j+2] = colors[i+2] / 255.0;

        vertsColsArray[j+3] = verts[i];
        vertsColsArray[j+4] = verts[i+1];
        vertsColsArray[j+5] = verts[i+2];
    }

    glGenBuffers( 2, BallsDataBuffers );

    glBindBuffer(GL_ARRAY_BUFFER, BallsDataBuffers[0]);
    glBufferData(GL_ARRAY_BUFFER, nVerts * 2 * sizeof(GLfloat), vertsColsArray, GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, BallsDataBuffers[1]);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, nTriangles * sizeof(GLuint), triangles, GL_STATIC_DRAW);

    triBallsToDraw = nTriangles;

    init = true;

    delete[] vertsColsArray;

    updateGL();
}

void MeshViewerWidget::loadLines(GLfloat* verts, GLfloat* colors, int nVerts, GLuint* lines, int nLines, QList<QPair<float, int> > es)
{
    GLfloat* linesColsArray = new GLfloat[nVerts * 2];

    for(int i = 0; i < nVerts; i = i + 3)
    {
        int j = 2 * i;
        linesColsArray[j] = colors[i] / 255.0;
        linesColsArray[j+1] = colors[i+1] / 255.0;
        linesColsArray[j+2] = colors[i+2] / 255.0;

        linesColsArray[j+3] = verts[i];
        linesColsArray[j+4] = verts[i+1];
        linesColsArray[j+5] = verts[i+2];
    }

    glGenBuffers( 2, LinesDataBuffers );


    glBindBuffer(GL_ARRAY_BUFFER, LinesDataBuffers[0]);
    glBufferData(GL_ARRAY_BUFFER, nVerts * 2 * sizeof(GLfloat), linesColsArray, GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, LinesDataBuffers[1]);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, nLines * sizeof(GLuint), lines, GL_STATIC_DRAW);

    linesToDraw = nLines;

    edgeSizes = es;

    delete[] linesColsArray;

    updateGL();
}

void MeshViewerWidget::loadPoints(GLfloat* verts, GLfloat* colors, int nVerts, GLuint* points, int nPoints, QList<QPair<float, int> > vs)
{
    GLfloat* pointsColsArray = new GLfloat[nVerts * 2];

    for(int i = 0; i < nVerts; i = i + 3)
    {
        int j = 2 * i;
        pointsColsArray[j] = colors[i] / 255.0;
        pointsColsArray[j+1] = colors[i+1] / 255.0;
        pointsColsArray[j+2] = colors[i+2] / 255.0;

        pointsColsArray[j+3] = verts[i];
        pointsColsArray[j+4] = verts[i+1];
        pointsColsArray[j+5] = verts[i+2];
    }

    glGenBuffers( 2, PointsDataBuffers );


    glBindBuffer(GL_ARRAY_BUFFER, PointsDataBuffers[0]);
    glBufferData(GL_ARRAY_BUFFER, nVerts * 2 * sizeof(GLfloat), pointsColsArray, GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, PointsDataBuffers[1]);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, nPoints * sizeof(GLuint), points, GL_STATIC_DRAW);

    pointsToDraw = nPoints;

    vertsSizes = vs;

    delete[] pointsColsArray;

    updateGL();
}


void MeshViewerWidget::paintGL()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode( GL_PROJECTION );
    glLoadMatrixd( projection_matrix_ );
    glMatrixMode( GL_MODELVIEW );
    glLoadMatrixd( modelview_matrix_ );

    if(triToDraw != 0)
    {
        glPolygonOffset(1.0, 2);

        // on charge le buffer 0 : une liste de vertex [r, g, b, x, y, z] (6 float)
        glBindBuffer(GL_ARRAY_BUFFER, TriDataBuffers[0]);

        // on charge la partie [r, g, b]
        glColorPointer( 3, GL_FLOAT, 6 * sizeof(float), 0 );
        // on charge la partie [x, y, z]
        glVertexPointer( 3, GL_FLOAT, 6 * sizeof(float), ((float*)NULL + (3)) );


        // on charge le buffer 1 : une liste d'ID [v0, v1, v2] (3 int)
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, TriDataBuffers[1]);

        glEnableClientState( GL_VERTEX_ARRAY );

        // affiche les faces en GL_COLOR_ARRAY
        glEnableClientState( GL_COLOR_ARRAY );
        glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );

        glEnable(GL_POLYGON_OFFSET_FILL);
        glDrawElements(GL_TRIANGLES, triToDraw, GL_UNSIGNED_INT, 0);
        glDisable(GL_POLYGON_OFFSET_FILL);

        glDisableClientState( GL_COLOR_ARRAY );
        glDisableClientState( GL_VERTEX_ARRAY );
    }

    // draw the balls
    if(triBallsToDraw != 0)
    {
        glPolygonOffset(1.0, 2);

        // on charge le buffer 0 : une liste de vertex [r, g, b, x, y, z] (6 float)
        glBindBuffer(GL_ARRAY_BUFFER, BallsDataBuffers[0]);

        // on charge la partie [r, g, b]
        glColorPointer( 3, GL_FLOAT, 6 * sizeof(float), 0 );
        // on charge la partie [x, y, z]
        glVertexPointer( 3, GL_FLOAT, 6 * sizeof(float), ((float*)NULL + (3)) );


        // on charge le buffer 1 : une liste d'ID [v0, v1, v2] (3 int)
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, BallsDataBuffers[1]);

        glEnableClientState( GL_VERTEX_ARRAY );

        // affiche les faces en GL_COLOR_ARRAY
        glEnableClientState( GL_COLOR_ARRAY );
        glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );

        glEnable(GL_POLYGON_OFFSET_FILL);
        glDrawElements(GL_TRIANGLES, triBallsToDraw, GL_UNSIGNED_INT, 0);
        glDisable(GL_POLYGON_OFFSET_FILL);

        glDisableClientState( GL_COLOR_ARRAY );
        glDisableClientState( GL_VERTEX_ARRAY );
    }


    if (linesToDraw != 0)
    {
        // on charge le buffer 0 : une liste de vertex [r, g, b, x, y, z] (6 float)
        glBindBuffer(GL_ARRAY_BUFFER, LinesDataBuffers[0]);

        // on charge la partie [r, g, b]
        glColorPointer( 3, GL_FLOAT, 6 * sizeof(float), 0 );
        // on charge la partie [x, y, z]
        glVertexPointer( 3, GL_FLOAT, 6 * sizeof(float), ((float*)NULL + (3)) );

        // on charge le buffer 1 : une liste d'ID [v0, v1] (2 int)
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, LinesDataBuffers[1]);

        glEnableClientState( GL_VERTEX_ARRAY );
        glEnableClientState( GL_COLOR_ARRAY );

        int cur = 0;
        for(int i = 0; i < edgeSizes.count(); i++)
        {
            glLineWidth(edgeSizes.at(i).first);
            glDrawElements(GL_LINES, edgeSizes.at(i).second*2, GL_UNSIGNED_INT, (GLvoid*)(sizeof(GLfloat) * cur));
            cur = cur + edgeSizes.at(i).second*2;
        }

        glDisableClientState( GL_COLOR_ARRAY );
        glDisableClientState( GL_VERTEX_ARRAY );
    }

    if (pointsToDraw != 0)
    {
        // on charge le buffer 0 : une liste de vertex [r, g, b, x, y, z] (6 float)
        glBindBuffer(GL_ARRAY_BUFFER, PointsDataBuffers[0]);

        // on charge la partie [r, g, b]
        glColorPointer( 3, GL_FLOAT, 6 * sizeof(float), 0 );
        // on charge la partie [x, y, z]
        glVertexPointer( 3, GL_FLOAT, 6 * sizeof(float), ((float*)NULL + (3)) );

        // on charge le buffer 1 : une liste d'ID [v0] (1 int)
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, PointsDataBuffers[1]);

        glEnableClientState( GL_VERTEX_ARRAY );
        glEnableClientState( GL_COLOR_ARRAY );

        int cur = 0;
        for(int i = 0; i < vertsSizes.count(); i++)
        {
            glPointSize(vertsSizes.at(i).first);
            glDrawElements(GL_POINTS, vertsSizes.at(i).second, GL_UNSIGNED_INT, (GLvoid*)(sizeof(GLfloat) * cur));
            cur = cur + vertsSizes.at(i).second;
        }

        glDisableClientState( GL_COLOR_ARRAY );
        glDisableClientState( GL_VERTEX_ARRAY );
    }
}


void MeshViewerWidget::set_scene_pos( const OpenMesh::Vec3f& _cog, float _radius )
{
    center_ = _cog;
    radius_ = _radius;
    update_projection_matrix();
    view_all();
}

void MeshViewerWidget::update_projection_matrix()
{
    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    gluPerspective(45.0, (GLfloat) width() / (GLfloat) height(), 0.01*radius_, 100.0*radius_);
//    glOrtho(0., (GLfloat)width(), 0., (GLfloat)height(), 0., 1.); // AldoGL
    glGetDoublev( GL_PROJECTION_MATRIX, projection_matrix_);
    glMatrixMode( GL_MODELVIEW );
}

void MeshViewerWidget::view_all()
{
    const GLdouble x = modelview_matrix_[0]*center_[0]  // M[0,0]*C[0]
            + modelview_matrix_[4]*center_[1]           // M[1,0]*C[1]
            + modelview_matrix_[8]*center_[2]           // M[2,0]*C[2]
            + modelview_matrix_[12];                    // M[3,0]*1

    const GLdouble y = modelview_matrix_[1]*center_[0]  // M[0,1]*C[0]
            + modelview_matrix_[5]*center_[1]           // M[1,1]*C[1]
            + modelview_matrix_[9]*center_[2]           // M[2,1]*C[2]
            + modelview_matrix_[13];                    // M[3,1]*1

    const GLdouble z = modelview_matrix_[2]*center_[0]  // M[0,2]*C[0]
            + modelview_matrix_[6]*center_[1]           // M[1,2]*C[1]
            + modelview_matrix_[10]*center_[2]          // M[2,2]*C[3]
            + modelview_matrix_[14]                     // M[3,2]*1
            + 3.0*radius_;

    translate(Vec3f(-x, -y, -z));
}

void MeshViewerWidget::mousePressEvent( QMouseEvent* _event )
{
    last_point_ok_ = map_to_sphere( last_point_2D_=_event->pos(), last_point_3D_ );
}

void MeshViewerWidget::mouseMoveEvent( QMouseEvent* _event )
{
    QPoint newPoint2D = _event->pos();

    Vec3f  newPoint3D;
    bool   newPoint_hitSphere = map_to_sphere( newPoint2D, newPoint3D );

    const float dx = newPoint2D.x() - last_point_2D_.x();
    const float dy = newPoint2D.y() - last_point_2D_.y();

    const float w  = width();
    const float h  = height();

    // zoom in/out with Ctrl+LeftButton
    if ( (_event->buttons() == (Qt::LeftButton+Qt::MidButton)) ||  (_event->buttons() == Qt::LeftButton && _event->modifiers() == Qt::ControlModifier))
    {
        const float value_y = radius_ * dy * 3.0 / h;
        translate(Vec3f(0.0, 0.0, value_y));
        // qDebug() << "zoom" << value_y;
    }
    // translate with MidButton or Alt+LeftButton
    else if ( (_event->buttons() == Qt::MidButton) || (_event->buttons() == Qt::LeftButton && _event->modifiers() == Qt::AltModifier) )
    {
        float z = modelview_matrix_[ 2]*center_[0]
                + modelview_matrix_[ 6]*center_[1]
                + modelview_matrix_[10]*center_[2]
                + modelview_matrix_[14];
        const float t = modelview_matrix_[ 3]*center_[0]
                + modelview_matrix_[ 7]*center_[1]
                + modelview_matrix_[11]*center_[2]
                + modelview_matrix_[15];
        z = -z / t;
        const float aspect     = w / h;
        const float near_plane = 0.01 * radius_;
        const float top        = tan(45.0f/2.0f*M_PI/180.0f) * near_plane;
        const float right      = aspect*top;
        translate(Vec3f( 2.0*dx/w*right/near_plane*z, -2.0*dy/h*top/near_plane*z, 0.0f));
//        qDebug() << "translate";

    }
    // rotate with Qt::LeftButton;
    else if (_event->buttons() == Qt::LeftButton)
    {
        if (last_point_ok_)
        {
            if ((newPoint_hitSphere = map_to_sphere(newPoint2D, newPoint3D)))
            {
                Vec3f axis = last_point_3D_ % newPoint3D;
                if (axis.sqrnorm() < 1e-7)
                    axis = Vec3f(1, 0, 0);
                else
                    axis.normalize();
                Vec3f d = last_point_3D_ - newPoint3D;
                float t = 0.5 * d.norm() / TRACKBALL_RADIUS;
                if (t < -1.0)
                    t = -1.0;
                else if (t > 1.0)
                    t = 1.0;
                float phi = 2.0 * asin(t);
                float angle = phi * 180.0 / M_PI;
                rotate(axis, angle);
            }
        }

    }
    last_point_2D_ = newPoint2D;
    last_point_3D_ = newPoint3D;
    last_point_ok_ = newPoint_hitSphere;
    updateGL();
}

void MeshViewerWidget::mouseReleaseEvent( QMouseEvent* /* _event */ )
{
    last_point_ok_ = false;
}

void MeshViewerWidget::wheelEvent(QWheelEvent* _event)
{
    float d = -(float)_event->delta() / 120.0 * 0.2 * radius_;
    translate(Vec3f(0.0, 0.0, d));
    updateGL();
    _event->accept();
}

bool MeshViewerWidget::map_to_sphere( const QPoint& _v2D, OpenMesh::Vec3f& _v3D )
{
    const double x =  (2.0*_v2D.x() - width())/width();
    const double y = -(2.0*_v2D.y() - height())/height();
    const double x2y2 = x*x + y*y;

    const double rsqr = TRACKBALL_RADIUS*TRACKBALL_RADIUS;
    _v3D[0] = x;
    _v3D[1] = y;
    if (x2y2 < 0.5*rsqr)
        _v3D[2] = sqrt(rsqr - x2y2);
    else
        _v3D[2] = 0.5*rsqr/sqrt(x2y2);

    return true;
}

void MeshViewerWidget::rotate( const OpenMesh::Vec3f& _axis, float _angle )
{
    const GLdouble x = modelview_matrix_[0]*center_[0] + modelview_matrix_[4]*center_[1] + modelview_matrix_[ 8]*center_[2] + modelview_matrix_[12];
    const GLdouble y = modelview_matrix_[1]*center_[0] + modelview_matrix_[5]*center_[1] + modelview_matrix_[ 9]*center_[2] + modelview_matrix_[13];
    const GLdouble z = modelview_matrix_[2]*center_[0] + modelview_matrix_[6]*center_[1] + modelview_matrix_[10]*center_[2] + modelview_matrix_[14];
    const Vec3f t(x, y, z);
    glLoadIdentity();
    glTranslatef(t[0], t[1], t[2]);
    glRotated( _angle, _axis[0], _axis[1], _axis[2]);
    glTranslatef(-t[0], -t[1], -t[2]);
    glMultMatrixd(modelview_matrix_);
    glGetDoublev(GL_MODELVIEW_MATRIX, modelview_matrix_);

    // debugging
//    qDebug() << "modelview_matrix_:";
//    for (int i = 0; i < 4; i++)
//        qDebug() << modelview_matrix_[4*i+0] << modelview_matrix_[4*i+1] << modelview_matrix_[4*i+2] << modelview_matrix_[4*i+3];
//    qDebug() << "projection_matrix_:";
//    for (int i = 0; i < 4; i++)
//        qDebug() << projection_matrix_[4*i+0] << projection_matrix_[4*i+1] << projection_matrix_[4*i+2] << projection_matrix_[4*i+3];
}


void MeshViewerWidget::setPOV(float x, float y)
{
    assert(0 <= x && x <= 2);
    assert(0 <= y && y <= 1);
    const double phi   = x * M_PI;
    const double theta = y * M_PI;

    modelview_matrix_[ 0] =  cosf(theta)*cosf(phi);
    modelview_matrix_[ 4] =  cosf(theta)*sinf(phi);
    modelview_matrix_[ 8] = -sinf(theta);
    modelview_matrix_[12] = 0;

    modelview_matrix_[ 1] = -sinf(phi);
    modelview_matrix_[ 5] =  cosf(phi);
    modelview_matrix_[ 9] = 0;
    modelview_matrix_[13] = 0;

    modelview_matrix_[ 2] = sinf(theta)*cosf(phi);
    modelview_matrix_[ 6] = sinf(theta)*sinf(phi);
    modelview_matrix_[10] = cosf(theta);
//    modelview_matrix_[14] = -3;

    modelview_matrix_[ 3] = 0;
    modelview_matrix_[ 7] = 0;
    modelview_matrix_[11] = 0;
    modelview_matrix_[15] = 1;
    updateGL();
}


void MeshViewerWidget::getPOV(double &x, double &y) const
{
    const float theta = acosf(modelview_matrix_[10]);
    const float phi = atan2f(modelview_matrix_[6], modelview_matrix_[2]);
    x = phi / M_PI;
    if (x < 0) x += 2.0;
    y = theta / M_PI;

    assert(0 <= x && x <= 2);
    assert(0 <= y && y <= 1);

//    qDebug() << "x, y =" << x << y;
}
