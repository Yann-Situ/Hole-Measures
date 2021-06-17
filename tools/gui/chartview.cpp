#include "chartview.h"

#include <QtCharts/QScatterSeries>
#include <QtCharts/QLegendMarker>
#include <QtGui/QImage>
#include <QtGui/QPainter>
#include <QtCore/QtMath>
#include <QDebug>

ChartView::ChartView(QWidget *parent) :
    QChartView(new QChart(), parent),
    m_pairs_0(0),
    m_pairs_1(0),
    m_pairs_2(0)
{

    setRenderHint(QPainter::Antialiasing);

//    chart()->setTitle("Thickness-breadth diagram");

    m_pairs_0 = new QScatterSeries();
    m_pairs_0->setName("dim 0");
    m_pairs_0->setMarkerSize(10.0);
    m_pairs_1 = new QScatterSeries();
    m_pairs_1->setName("dim 1");
    m_pairs_1->setMarkerSize(10.0);
    m_pairs_2 = new QScatterSeries();
    m_pairs_2->setName("dim 2");
    m_pairs_2->setMarkerSize(10.0);

    chart()->addSeries(m_pairs_0);
    chart()->addSeries(m_pairs_1);
    chart()->addSeries(m_pairs_2);
    chart()->createDefaultAxes();
    chart()->legend()->setMarkerShape(QLegend::MarkerShapeFromSeries);
    chart()->margins().setLeft(5);
    chart()->margins().setBottom(5);
    chart()->margins().setRight(5);
    chart()->margins().setTop(5);

    connect(m_pairs_0, &QScatterSeries::clicked, this, &ChartView::handleClickedPoint);
    connect(m_pairs_1, &QScatterSeries::clicked, this, &ChartView::handleClickedPoint);
    connect(m_pairs_2, &QScatterSeries::clicked, this, &ChartView::handleClickedPoint);
}


void ChartView::load_pairs(std::vector<int> dim, std::vector<double> x, std::vector<double> y)
{
    // reset the diagram
    m_pairs_0->clear();
    m_pairs_1->clear();
    m_pairs_2->clear();

    // add pairs to the series
    double max_value = 0;
    for (std::size_t i = 0; i < dim.size(); ++i)
    {
        max_value = std::max(max_value, x.at(i));
        max_value = std::max(max_value, y.at(i));
        if (dim.at(i) == 0)
            m_pairs_0->append(x.at(i), y.at(i));
        else if (dim.at(i) == 1)
            m_pairs_1->append(x.at(i), y.at(i));
        else if (dim.at(i) == 2)
            m_pairs_2->append(x.at(i), y.at(i));
    }
    chart()->createDefaultAxes();
    chart()->axes(Qt::Horizontal).first()->setRange(0, 1.1 * max_value);
    chart()->axes(Qt::Vertical  ).first()->setRange(0, 1.1 * max_value);
}


void ChartView::keyPressEvent(QKeyEvent *event)
{
    switch (event->key()) {
    case Qt::Key_Plus:
        chart()->zoom(1.5);
        break;
    case Qt::Key_Minus:
        chart()->zoom(0.7);
        break;
    case Qt::Key_Left:
        chart()->scroll(-10, 0);
        break;
    case Qt::Key_Right:
        chart()->scroll(10, 0);
        break;
    case Qt::Key_Up:
        chart()->scroll(0, 10);
        break;
    case Qt::Key_Down:
        chart()->scroll(0, -10);
        break;
    default:
        QGraphicsView::keyPressEvent(event);
        break;
    }
//    qDebug() << chart()->plotArea();
}

void ChartView::handleClickedPoint(const QPointF &point)
{
    QPointF clickedPoint = point;
    // Find the closest point from series 1
    int dim_closest = 0;
    int pos_closest = 0;
    qreal distance(INT_MAX);
    // dimension 0
    auto points = m_pairs_0->pointsVector();
    for (int i = 0; i < points.size(); ++i)
    {
        qreal currentDistance = dist(points[i], clickedPoint);
        if (currentDistance < distance) {
            distance = currentDistance;
            dim_closest = 0;
            pos_closest = i;
        }
    }
    // dimension 1
    points = m_pairs_1->pointsVector();
    for (int i = 0; i < points.size(); ++i)
    {
        qreal currentDistance = dist(points[i], clickedPoint);
        if (currentDistance < distance) {
            distance = currentDistance;
            dim_closest = 1;
            pos_closest = i;
        }
    }
    // dimension 2
    points = m_pairs_2->pointsVector();
    for (int i = 0; i < points.size(); ++i)
    {
        qreal currentDistance = dist(points[i], clickedPoint);
        if (currentDistance < distance) {
            distance = currentDistance;
            dim_closest = 2;
            pos_closest = i;
        }
    }
    qDebug() << "selected point: " << dim_closest << " " << pos_closest;
    emit pairSelected(dim_closest, pos_closest);

    // Remove the closes point from series 1 and append it to series 2
//    m_scatter->remove(closest);
//    m_scatter2->append(closest);
}


qreal ChartView::dist(const QPointF &p1, const QPointF &p2) const
{
    return qSqrt((p1.x() - p2.x()) * (p1.x() - p2.x())
                 + (p1.y() - p2.y()) * (p1.y() - p2.y()));
}
