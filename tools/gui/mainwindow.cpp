#include "mainwindow.h"
#include "ui_mainwindow.h"

#include <QMessageBox>
//#include <QtCharts/QScatterSeries>

#include <fstream>	// ifstream

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    m_thickness(false),
    max_persistence(0.0),
    current_min_persistence(0.0),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);

    QObject::connect(ui->actionOpenMesh, SIGNAL(triggered()), this, SLOT(read_mesh()));
    QObject::connect(ui->actionOpenTB  , SIGNAL(triggered()), this, SLOT(read_tb()  ));
    QObject::connect(ui->chartView_diagram, &ChartView::pairSelected, this, &MainWindow::selectPair);
}

MainWindow::~MainWindow()
{
    delete ui;
}


void MainWindow::on_pushButton_test_clicked()
{
//    TBball ball;
//    ball.dim = 1;
//    ball.t_radius = 1;
//    ball.t_x = 0; ball.t_y = 0; ball.t_z = 0;
//    show_tb_ball(ball);

    print_tb_pairs();
}


void MainWindow::on_checkBox_wireframe_stateChanged()
{
    displayMesh(&mesh);
}

void MainWindow::on_comboBox_TBchoice_currentIndexChanged()
{
    display_balls();
}

void MainWindow::on_checkBox_allballs_stateChanged()
{
    display_balls();
}

void MainWindow::on_spinBox_hole_valueChanged(int i)
{
    ui->checkBox_allballs->setCheckState(Qt::Unchecked);
    int dim = m_tb_balls.at(i).dim;
    if (dim==0)
        ui->checkBox_dim0->setCheckState(Qt::Checked);
    else if (dim==1)
        ui->checkBox_dim1->setCheckState(Qt::Checked);
    else if (dim==2)
        ui->checkBox_dim2->setCheckState(Qt::Checked);
    ui->spinbox_min_persistence->setValue(m_tb_balls.at(i).persistence*0.999);
    display_balls();
}

void MainWindow::on_checkBox_dim0_stateChanged()
{
    display_balls();
}

void MainWindow::on_checkBox_dim1_stateChanged()
{
    display_balls();
}

void MainWindow::on_checkBox_dim2_stateChanged()
{
    display_balls();
}

void MainWindow::on_slider_min_persistence_valueChanged()
{
    const int v = ui->slider_min_persistence->value();
    const int mini = ui->slider_min_persistence->minimum();
    const int maxi = ui->slider_min_persistence->maximum();
    const float range = (maxi-mini);
    if (range <= 0.0)
        current_min_persistence = 0.0;
    else
        current_min_persistence = (max_persistence*(v-mini))/(maxi-mini);
    ui->spinbox_min_persistence->blockSignals(True);
    ui->spinbox_min_persistence->setValue(current_min_persistence);
    ui->spinbox_min_persistence->blockSignals(False);
    display_balls();
}

void MainWindow::on_spinbox_min_persistence_valueChanged()
{
    const int mini = ui->slider_min_persistence->minimum();
    const int maxi = ui->slider_min_persistence->maximum();
    const double v = ui->spinbox_min_persistence->value();
    if (v <= 0.0)
        current_min_persistence = 0.0;
    else if (v >= max_persistence)
        current_min_persistence = max_persistence;
    else
        current_min_persistence = v;

    ui->slider_min_persistence->blockSignals(True);
    ui->slider_min_persistence->setValue(int(1.0*current_min_persistence/max_persistence * (maxi-mini) + mini));
    ui->slider_min_persistence->blockSignals(False);
    display_balls();
}

/**
 * @brief MainWindow::read_mesh
 * Load a mesh file and display it
 */
void MainWindow::read_mesh()
{
    // fenêtre de sélection des fichiers
    const QString fileName = QFileDialog::getOpenFileName(this, tr("Open Mesh"), "../../data", tr("Mesh Files (*.obj, *.off)"));
    m_filebasename = QFileInfo(fileName).baseName();
    setWindowTitle("Hole Measures GUI - " + QFileInfo(fileName).fileName());

    // if no file was chosen, stop here
    if (fileName.isEmpty()) return;

    // chargement du fichier .obj dans la variable globale "mesh"
    OpenMesh::IO::read_mesh(mesh, fileName.toUtf8().constData());

    mesh.update_normals();

    // initialisation des couleurs et épaisseurs (sommets et arêtes) du mesh
    resetAllColorsAndThickness(&mesh);

    // on affiche le maillage
    displayMesh(&mesh);
    ui->checkBox_wireframe->setEnabled(true);
    qDebug() << fileName << "loaded";
}

/**
 * @brief MainWindow::read_tb
 * Read a TB file containing the dimension, and TB balls of each hole
 */
void MainWindow::read_tb()
{
    const QString fileName = QFileDialog::getOpenFileName(this, tr("Open TB file"), "../../data", tr("TB files (*.tb)"));
    m_filebasename = QFileInfo(fileName).baseName();

    // if no file was chosen, stop here
    if (fileName.isEmpty()) return;

    std::ifstream file;
    file.open(fileName.toStdString().c_str(), std::ios::in);
    if (!file.good())
    {
        QMessageBox::warning(this, "File", "Cannot open the TB file");
        exit(EXIT_FAILURE);
    }
    std::string line;
    m_tb_balls.clear();
    max_persistence = 0.0;
    while (getline(file, line))
    {
        TBball ball;
        std::istringstream iss(line);
        qDebug() << QString::fromStdString(line);
        iss >> ball.dim;
        iss >> ball.t_radius;
        iss >> ball.t_x >> ball.t_y >> ball.t_z;
        iss >> ball.b_radius;
        iss >> ball.b_x >> ball.b_y >> ball.b_z;
        ball.persistence = ball.t_radius + ball.b_radius;
        m_tb_balls.push_back(ball);
        max_persistence = (ball.persistence > max_persistence) ? ball.persistence : max_persistence;
        qDebug() << ball.t_radius << " " << ball.b_radius;
    }
    file.close();

    m_displayed_balls.resize(m_tb_balls.size(), false);

    qDebug() << fileName << "loaded";

    ui->comboBox_TBchoice->setEnabled(true);
    ui->spinBox_hole->setMaximum(m_tb_balls.size()-1);
    ui->spinBox_hole->setEnabled(true);
    ui->checkBox_allballs->setEnabled(true);
    ui->checkBox_dim0->setEnabled(true);
    ui->checkBox_dim1->setEnabled(true);
    ui->checkBox_dim2->setEnabled(true);

    ui->slider_min_persistence->setEnabled(true);
    ui->spinbox_min_persistence->setEnabled(true);
    ui->spinbox_min_persistence->setValue(current_min_persistence);
    ui->spinbox_min_persistence->setMinimum(0.0);
    ui->spinbox_min_persistence->setMaximum(max_persistence);
    ui->spinbox_min_persistence->setSingleStep(max_persistence/100.0);

    display_balls();
    print_tb_pairs();
}

void MainWindow::show_tb_balls(const std::vector<TBball>& balls)
{
    const int divs = 16;
    const double inc = 3.14159265359 / (divs);
    // count faces to display

    double b_r, b_x, b_y, b_z;
    std::vector<int> color(3, 0);

    // precompute the vertices coordinates
    std::vector<std::vector<std::vector<double>>> vertices_vec;
    std::vector<bool> thickness_balls;
    std::vector<size_t> ball_hole_index;
    std::vector<double> v(3, 0);

    ui->text_HoleInfo->clear();
    for (size_t k = 0; k < balls.size(); k++) {
        const TBball& ball = balls.at(k);
        if ((ball.dim==0 && !ui->checkBox_dim0->isChecked()) ||
            (ball.dim==1 && !ui->checkBox_dim1->isChecked()) ||
            (ball.dim==2 && !ui->checkBox_dim2->isChecked()) )
            continue;

        if (ball.persistence < current_min_persistence)
            continue;

        if (ui->comboBox_TBchoice->currentIndex() == 0 || ui->comboBox_TBchoice->currentIndex() == 2)
        {
            b_r = ball.t_radius;
            b_x = ball.t_x;
            b_y = ball.t_y;
            b_z = ball.t_z;

            std::vector<std::vector<double>> vertices;
            v[0] = b_x; v[1] = b_y; v[2] = b_z+b_r; vertices.push_back(v);
            v[0] = b_x; v[1] = b_y; v[2] = b_z-b_r; vertices.push_back(v);
            for (int t = 1; t < divs; t++)	// t in {pi/4, 2*pi/4, 3*pi/4}
            {
                const double r_sin_t = b_r * sin(t * inc);
                const double z = b_r * cos(t * inc);
                for(int p = 0; p < 2*divs; p++)	// p in {0, pi/4, ..., 7*pi/4}
                {
                    const double x = r_sin_t * cos(p * inc);
                    const double y = r_sin_t * sin(p * inc);
                    v[0] = b_x+x; v[1] = b_y+y; v[2] = b_z+z; vertices.push_back(v);
                }
            }
            vertices_vec.push_back(vertices);
            thickness_balls.push_back(true);
            ball_hole_index.push_back(k);
        }
        if (ui->comboBox_TBchoice->currentIndex() == 1 || ui->comboBox_TBchoice->currentIndex() == 2)
        {
            b_r = ball.b_radius;
            b_x = ball.b_x;
            b_y = ball.b_y;
            b_z = ball.b_z;

            std::vector<std::vector<double>> vertices;
            v[0] = b_x; v[1] = b_y; v[2] = b_z+b_r; vertices.push_back(v);
            v[0] = b_x; v[1] = b_y; v[2] = b_z-b_r; vertices.push_back(v);
            for (int t = 1; t < divs; t++)	// t in {pi/4, 2*pi/4, 3*pi/4}
            {
                const double r_sin_t = b_r * sin(t * inc);
                const double z = b_r * cos(t * inc);
                for(int p = 0; p < 2*divs; p++)	// p in {0, pi/4, ..., 7*pi/4}
                {
                    const double x = r_sin_t * cos(p * inc);
                    const double y = r_sin_t * sin(p * inc);
                    v[0] = b_x+x; v[1] = b_y+y; v[2] = b_z+z; vertices.push_back(v);
                }
            }
            vertices_vec.push_back(vertices);
            thickness_balls.push_back(false);
            ball_hole_index.push_back(k);
        }

        ui->text_HoleInfo->append("dim: " +
                QString::number(ball.dim) + " (t,b): (" +
                QString::number(ball.t_radius) + ", " + QString::number(ball.b_radius) + ")"
        );
    }

    /* Load faces */
    unsigned displayed_faces = vertices_vec.size()*(2*2*divs + 2*2*divs*(divs-2));
    GLuint* triIndiceArray = new GLuint[displayed_faces * 3];
    GLfloat* triCols = new GLfloat[displayed_faces * 3 * 3];  // color of each vertex of each face
    GLfloat* triVerts = new GLfloat[displayed_faces * 3 * 3]; // coordinates of each vertex of each face
    int i = 0;

    int k = 0;
    for (std::vector<std::vector<double>> vertices : vertices_vec)
    {
        if (thickness_balls.at(k))
        {
            if (balls.at(ball_hole_index[k]).t_radius > 0.0)
            {
                color[0] = 250;
                color[1] = 40;
                color[2] = 80;
            }
            else
            {
                color[0] = 250;
                color[1] = 100;
                color[2] = 200;
            }
        }
        else
        {
            if (balls.at(ball_hole_index[k]).b_radius > 0.0)
            {
                color[0] = 40;
                color[1] = 80;
                color[2] = 250;
            }
            else
            {
                color[0] = 100;
                color[1] = 200;
                color[2] = 250;
            }
        }
        // Put triangular faces
        for(int p = 0; p < 2*divs; p++)
        {
            const int a = 2 + p;
            const int b = 2 + (p+1)%(2*divs);

            triCols[3*i+0] = color[0]; // RGB components of its color
            triCols[3*i+1] = color[1];
            triCols[3*i+2] = color[2];
            triVerts[3*i+0] = vertices.at(0)[0];
            triVerts[3*i+1] = vertices.at(0)[1];
            triVerts[3*i+2] = vertices.at(0)[2];
            triIndiceArray[i] = i;
            i++;
            triCols[3*i+0] = color[0];
            triCols[3*i+1] = color[1];
            triCols[3*i+2] = color[2];
            triVerts[3*i+0] = vertices.at(a)[0];
            triVerts[3*i+1] = vertices.at(a)[1];
            triVerts[3*i+2] = vertices.at(a)[2];
            triIndiceArray[i] = i;
            i++;
            triCols[3*i+0] = color[0];
            triCols[3*i+1] = color[1];
            triCols[3*i+2] = color[2];
            triVerts[3*i+0] = vertices.at(b)[0];
            triVerts[3*i+1] = vertices.at(b)[1];
            triVerts[3*i+2] = vertices.at(b)[2];
            triIndiceArray[i] = i;
            i++;

            const int c = vertices.size()-2*divs + (p+1)%(2*divs);
            const int d = vertices.size()-2*divs + p;
            triCols[3*i+0] = color[0];
            triCols[3*i+1] = color[1];
            triCols[3*i+2] = color[2];
            triVerts[3*i+0] = vertices.at(1)[0];
            triVerts[3*i+1] = vertices.at(1)[1];
            triVerts[3*i+2] = vertices.at(1)[2];
            triIndiceArray[i] = i;
            i++;
            triCols[3*i+0] = color[0];
            triCols[3*i+1] = color[1];
            triCols[3*i+2] = color[2];
            triVerts[3*i+0] = vertices.at(c)[0];
            triVerts[3*i+1] = vertices.at(c)[1];
            triVerts[3*i+2] = vertices.at(c)[2];
            triIndiceArray[i] = i;
            i++;
            triCols[3*i+0] = color[0];
            triCols[3*i+1] = color[1];
            triCols[3*i+2] = color[2];
            triVerts[3*i+0] = vertices.at(d)[0];
            triVerts[3*i+1] = vertices.at(d)[1];
            triVerts[3*i+2] = vertices.at(d)[2];
            triIndiceArray[i] = i;
            i++;
        }
        // put square faces
        for (int t = 1; t <= divs-2; t++)
        {
            for (int p = 0; p < 2*divs; p++)
            {
                const int a = 2 + 2*divs*(t-1) + p;
                const int b = 2 + 2*divs*t + p;
                const int c = 2 + 2*divs*t + (p+1)%(2*divs);
                const int d = 2 + 2*divs*(t-1) + (p+1)%(2*divs);
                triCols[3*i+0] = color[0];
                triCols[3*i+1] = color[1];
                triCols[3*i+2] = color[2];
                triVerts[3*i+0] = vertices.at(a)[0];
                triVerts[3*i+1] = vertices.at(a)[1];
                triVerts[3*i+2] = vertices.at(a)[2];
                triIndiceArray[i] = i;
                i++;
                triCols[3*i+0] = color[0];
                triCols[3*i+1] = color[1];
                triCols[3*i+2] = color[2];
                triVerts[3*i+0] = vertices.at(b)[0];
                triVerts[3*i+1] = vertices.at(b)[1];
                triVerts[3*i+2] = vertices.at(b)[2];
                triIndiceArray[i] = i;
                i++;
                triCols[3*i+0] = color[0];
                triCols[3*i+1] = color[1];
                triCols[3*i+2] = color[2];
                triVerts[3*i+0] = vertices.at(c)[0];
                triVerts[3*i+1] = vertices.at(c)[1];
                triVerts[3*i+2] = vertices.at(c)[2];
                triIndiceArray[i] = i;
                i++;

                triCols[3*i+0] = color[0];
                triCols[3*i+1] = color[1];
                triCols[3*i+2] = color[2];
                triVerts[3*i+0] = vertices.at(a)[0];
                triVerts[3*i+1] = vertices.at(a)[1];
                triVerts[3*i+2] = vertices.at(a)[2];
                triIndiceArray[i] = i;
                i++;
                triCols[3*i+0] = color[0];
                triCols[3*i+1] = color[1];
                triCols[3*i+2] = color[2];
                triVerts[3*i+0] = vertices.at(c)[0];
                triVerts[3*i+1] = vertices.at(c)[1];
                triVerts[3*i+2] = vertices.at(c)[2];
                triIndiceArray[i] = i;
                i++;
                triCols[3*i+0] = color[0];
                triCols[3*i+1] = color[1];
                triCols[3*i+2] = color[2];
                triVerts[3*i+0] = vertices.at(d)[0];
                triVerts[3*i+1] = vertices.at(d)[1];
                triVerts[3*i+2] = vertices.at(d)[2];
                triIndiceArray[i] = i;
                i++;
            }
        }

        k++;
    }

    ui->displayWidget->loadBalls(triVerts, triCols, displayed_faces * 3 * 3, triIndiceArray, displayed_faces * 3);
//    qDebug() << "One ball shown";
    /*ui->statusBar->showMessage(
                "dim: " + QString::number(ball.dim) +
                " (t,b): (" + QString::number(ball.t_radius) + ", " + QString::number(ball.b_radius) + ")");*/
}

void MainWindow::display_balls()
{
    if (ui->checkBox_allballs->isChecked())
    {
        show_tb_balls(m_tb_balls);
    }
    else
    {
        std::vector<TBball> balls;
        balls.push_back(m_tb_balls.at(ui->spinBox_hole->value()));
        show_tb_balls(balls);
    }
}

/* **** fonctions supplémentaires **** */
// permet d'initialiser les couleurs et les épaisseurs des élements du maillage
void MainWindow::resetAllColorsAndThickness(MyMesh* _mesh)
{
    // set thickness (1) and color (black) for vertices
    for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
    {
        _mesh->data(*curVert).thickness = 2;
        _mesh->set_color(*curVert, MyMesh::Color(0, 0, 0));
    }
    // set color (light gray) for faces
    for (MyMesh::FaceIter curFace = _mesh->faces_begin(); curFace != _mesh->faces_end(); curFace++)
    {
        _mesh->set_color(*curFace, MyMesh::Color(220, 220, 220));
    }
    // set thickness (1) and color (dark gray) for edges
    for (MyMesh::EdgeIter curEdge = _mesh->edges_begin(); curEdge != _mesh->edges_end(); curEdge++)
    {
        _mesh->data(*curEdge).thickness = 1;
        _mesh->set_color(*curEdge, MyMesh::Color(90, 90, 90));
    }
}

/**
 * @brief Affiche le maillage. Ecrit par Arnaud Polette.
 * @param _mesh
 * @param isTemperatureMap
 * @param mapRange
 */
void MainWindow::displayMesh(MyMesh* _mesh, bool isTemperatureMap, float mapRange) const
{
    // display the faces
    if (!ui->checkBox_wireframe->isChecked())
    {
        // count faces to display
        unsigned displayed_faces = 0;
        for (MyMesh::FaceIter curFace = _mesh->faces_begin(); curFace != _mesh->faces_end(); curFace++)
        {
            if (_mesh->color(*curFace) != MyMesh::Color(0,255,0)) /// @note What??
                displayed_faces++;
        }

        GLuint* triIndiceArray = new GLuint[displayed_faces * 3];
        GLfloat* triCols = new GLfloat[displayed_faces * 3 * 3];  // color of each vertex of each face
        GLfloat* triVerts = new GLfloat[displayed_faces * 3 * 3]; // coordinates of each vertex of each face

        int i = 0;

        if (isTemperatureMap) // if vertices have values (?)
        {
            QVector<float> values;

            if (mapRange == -1)
            {
                for (MyMesh::VertexIter curVert = _mesh->vertices_begin(); curVert != _mesh->vertices_end(); curVert++)
                    values.append(fabs(_mesh->data(*curVert).value));
                std::sort(values.begin(), values.end());
                //            qSort(values);
                mapRange = values.at(values.size()*0.8); // almost the larger value of the vertices
                qDebug() << "mapRange" << mapRange;
            }

            const float range = mapRange;
            MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
            MyMesh::ConstFaceVertexIter fvIt;

            for (; fIt != fEnd; ++fIt) // for each face
            {
                fvIt = _mesh->cfv_iter(*fIt);
                if (_mesh->data(*fvIt).value > 0)
                {
                    triCols[3*i+0] = 255;
                    triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);
                    triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);
                }
                else
                {
                    triCols[3*i+2] = 255;
                    triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);
                    triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);
                }
                triVerts[3*i+0] = _mesh->point(*fvIt)[0];
                triVerts[3*i+1] = _mesh->point(*fvIt)[1];
                triVerts[3*i+2] = _mesh->point(*fvIt)[2];
                triIndiceArray[i] = i;

                i++; ++fvIt;
                if(_mesh->data(*fvIt).value > 0)
                {
                    triCols[3*i+0] = 255;
                    triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);
                    triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);
                }
                else
                {
                    triCols[3*i+2] = 255;
                    triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);
                    triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);
                }
                triVerts[3*i+0] = _mesh->point(*fvIt)[0];
                triVerts[3*i+1] = _mesh->point(*fvIt)[1];
                triVerts[3*i+2] = _mesh->point(*fvIt)[2];
                triIndiceArray[i] = i;

                i++; ++fvIt;
                if(_mesh->data(*fvIt).value > 0)
                {
                    triCols[3*i+0] = 255;
                    triCols[3*i+1] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);
                    triCols[3*i+2] = 255 - std::min((_mesh->data(*fvIt).value/range) * 255.0, 255.0);
                }
                else
                {
                    triCols[3*i+2] = 255;
                    triCols[3*i+1] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);
                    triCols[3*i+0] = 255 - std::min((-_mesh->data(*fvIt).value/range) * 255.0, 255.0);
                }
                triVerts[3*i+0] = _mesh->point(*fvIt)[0];
                triVerts[3*i+1] = _mesh->point(*fvIt)[1];
                triVerts[3*i+2] = _mesh->point(*fvIt)[2];
                triIndiceArray[i] = i;

                i++;
            }
        }
        else
        {
            /* Load faces */
            MyMesh::ConstFaceIter fIt(_mesh->faces_begin()), fEnd(_mesh->faces_end());
            MyMesh::ConstFaceVertexIter fvIt;
            for (; fIt != fEnd; ++fIt) // for each face not green (?)
            {
                if (_mesh->color(*fIt) != MyMesh::Color(0,255,0))
                {
                    fvIt = _mesh->cfv_iter(*fIt); // one vertex of the face
                    triCols[3*i+0] = _mesh->color(*fIt)[0]; // RGB components of its color
                    triCols[3*i+1] = _mesh->color(*fIt)[1];
                    triCols[3*i+2] = _mesh->color(*fIt)[2];
                    triVerts[3*i+0] = _mesh->point(*fvIt)[0]; // coordinates of the first vertex
                    triVerts[3*i+1] = _mesh->point(*fvIt)[1];
                    triVerts[3*i+2] = _mesh->point(*fvIt)[2];
                    triIndiceArray[i] = i;

                    i++; ++fvIt; // idem for next vertex of the face
                    triCols[3*i+0] = _mesh->color(*fIt)[0];
                    triCols[3*i+1] = _mesh->color(*fIt)[1];
                    triCols[3*i+2] = _mesh->color(*fIt)[2];
                    triVerts[3*i+0] = _mesh->point(*fvIt)[0];
                    triVerts[3*i+1] = _mesh->point(*fvIt)[1];
                    triVerts[3*i+2] = _mesh->point(*fvIt)[2];
                    triIndiceArray[i] = i;

                    i++; ++fvIt;
                    triCols[3*i+0] = _mesh->color(*fIt)[0];
                    triCols[3*i+1] = _mesh->color(*fIt)[1];
                    triCols[3*i+2] = _mesh->color(*fIt)[2];
                    triVerts[3*i+0] = _mesh->point(*fvIt)[0];
                    triVerts[3*i+1] = _mesh->point(*fvIt)[1];
                    triVerts[3*i+2] = _mesh->point(*fvIt)[2];
                    triIndiceArray[i] = i;
                    i++;
                }
            }
        }

        ui->displayWidget->loadMesh(triVerts, triCols, displayed_faces * 3 * 3, triIndiceArray, displayed_faces * 3);

        delete[] triIndiceArray;
        delete[] triCols;
        delete[] triVerts;
    }
    else
    {
        GLuint* triIndiceArray = new GLuint[0 * 3];
        GLfloat* triCols = new GLfloat[0 * 3 * 3];  // color of each vertex of each face
        GLfloat* triVerts = new GLfloat[0 * 3 * 3]; // coordinates of each vertex of each face

        ui->displayWidget->loadMesh(triVerts, triCols, 0 * 3 * 3, triIndiceArray, 0 * 3);

        delete[] triIndiceArray;
        delete[] triCols;
        delete[] triVerts;
    }



    /* Load edges */
    GLuint* linesIndiceArray = new GLuint[_mesh->n_edges() * 2];
    GLfloat* linesCols = new GLfloat[_mesh->n_edges() * 2 * 3];
    GLfloat* linesVerts = new GLfloat[_mesh->n_edges() * 2 * 3];

    // put edges in a map grouped by thickness
    int i = 0;
    QHash<float, QList<int> > edgesIDbyThickness;
    for (MyMesh::EdgeIter eit = _mesh->edges_begin(); eit != _mesh->edges_end(); ++eit) // for each edge
    {
        float t = _mesh->data(*eit).thickness;
        if (t > 0)
        {
            if (!edgesIDbyThickness.contains(t))
                edgesIDbyThickness[t] = QList<int>();
            edgesIDbyThickness[t].append((*eit).idx());
        }
    }
    QHashIterator<float, QList<int> > it(edgesIDbyThickness);
    QList<QPair<float, int> > edgeSizes;
    while (it.hasNext()) // for each edge
    {
        it.next();

        for (int e = 0; e < it.value().size(); e++)
        {
            const int eidx = it.value().at(e);

            MyMesh::VertexHandle vh1 = _mesh->to_vertex_handle(_mesh->halfedge_handle(_mesh->edge_handle(eidx), 0));
            if (!_mesh->face_handle(_mesh->halfedge_handle(_mesh->edge_handle(eidx), 0)).is_valid() &&
                    !_mesh->face_handle(_mesh->halfedge_handle(_mesh->edge_handle(eidx), 1)).is_valid()) // NOTE ??
                continue;

            linesVerts[3*i+0] = _mesh->point(vh1)[0]; // coordinates of one vertex of the edge
            linesVerts[3*i+1] = _mesh->point(vh1)[1];
            linesVerts[3*i+2] = _mesh->point(vh1)[2];
            linesCols[3*i+0] = _mesh->color(_mesh->edge_handle(eidx))[0]; // RGB components of the edge
            linesCols[3*i+1] = _mesh->color(_mesh->edge_handle(eidx))[1];
            linesCols[3*i+2] = _mesh->color(_mesh->edge_handle(eidx))[2];
            linesIndiceArray[i] = i;
            i++;

            MyMesh::VertexHandle vh2 = _mesh->from_vertex_handle(_mesh->halfedge_handle(_mesh->edge_handle(eidx), 0)); // NOTE idem to vh1!
            linesVerts[3*i+0] = _mesh->point(vh2)[0];
            linesVerts[3*i+1] = _mesh->point(vh2)[1];
            linesVerts[3*i+2] = _mesh->point(vh2)[2];
            linesCols[3*i+0] = _mesh->color(_mesh->edge_handle(eidx))[0];
            linesCols[3*i+1] = _mesh->color(_mesh->edge_handle(eidx))[1];
            linesCols[3*i+2] = _mesh->color(_mesh->edge_handle(eidx))[2];
            linesIndiceArray[i] = i;
            i++;
        }
        edgeSizes.append(qMakePair(it.key(), it.value().size()));
    }

    ui->displayWidget->loadLines(linesVerts, linesCols, i * 3, linesIndiceArray, i, edgeSizes);

    delete[] linesIndiceArray;
    delete[] linesCols;
    delete[] linesVerts;
}


/**
 * @brief MainWindow::displayPath
 * I may use this for displaying the balls...
 */
//void MainWindow::displayPath(vector<CPoint3D> resultpoints)
//{
//    GLuint* linesIndiceArray = new GLuint[(resultpoints.size()+2) * 2];
//    GLfloat* linesCols = new GLfloat[(resultpoints.size()+2) * 2 * 3];
//    GLfloat* linesVerts = new GLfloat[(resultpoints.size()+2) * 2 * 3];

//    QList<QPair<float, int> > edgeSizes;
//    std::size_t i = 0;
//    int j = 0;
//    while (i < resultpoints.size()*2 - 2)
//    {
//        linesVerts[3*i+0] = resultpoints[j].x;
//        linesVerts[3*i+1] = resultpoints[j].y;
//        linesVerts[3*i+2] = resultpoints[j].z;
//        linesCols[3*i+0] = 0;
//        linesCols[3*i+1] = 100;
//        linesCols[3*i+2] = 0;
//        linesIndiceArray[i] = i;
//        edgeSizes.append(qMakePair(10, i));
//        i++;
//        j++;

//        linesVerts[3*i+0] = resultpoints[j].x;
//        linesVerts[3*i+1] = resultpoints[j].y;
//        linesVerts[3*i+2] = resultpoints[j].z;
//        linesCols[3*i+0] = 0;
//        linesCols[3*i+1] = 100;
//        linesCols[3*i+2] = 0;
//        linesIndiceArray[i] = i;
//        edgeSizes.append(qMakePair(10, i));
//        i++;
//    }

//    ui->displayWidget->loadPath(linesVerts, linesCols, i * 3, linesIndiceArray, i, edgeSizes);

//    delete[] linesIndiceArray;
//    delete[] linesCols;
//    delete[] linesVerts;
//}

/**
 * @brief MainWindow::print_tb_pairs
 * Use the TBballs of the TB file to update the TB diagram
 */
void MainWindow::print_tb_pairs()
{
    std::vector<int> dim;
    std::vector<double> x, y;
    for (TBball tb_ball : m_tb_balls)
    {
        dim.push_back(tb_ball.dim);
        x.push_back(tb_ball.t_radius);
        y.push_back(tb_ball.b_radius);
    }
    ui->chartView_diagram->load_pairs(dim, x, y);
}

void MainWindow::selectPair(int dim, int pos)
{
    for (std::size_t i = 0; i < m_tb_balls.size(); ++i)
    {
        if (m_tb_balls.at(i).dim == dim)
        {
            if (pos > 0)
                pos--;
            else
            {
                ui->spinBox_hole->setValue(i);
                return;
            }
        }
    }

}
