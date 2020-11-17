#include "ToolsConfiguration.h"

ToolsConfiguration::ToolsConfiguration(QWidget *parent) :
	QDialog(parent), ui(new Ui::ToolsConfig)
{
	ui->setupUi(this);

	InitializeTools();
	InitializeTable();

	connect(ui->listWidget, SIGNAL(itemClicked(QListWidgetItem*)), this, SLOT(onChangedToolSelection(QListWidgetItem*)));
	connect(ui->cancelButton, SIGNAL(clicked()), parent, SLOT(close()));
	connect(ui->addButton, SIGNAL(clicked()), this, SLOT(onActionAddTool()));
}


ToolsConfiguration::~ToolsConfiguration()
{
	delete ui;
}

void ToolsConfiguration::InitializeTools()
{
	currentListItem = nullptr;
	QListWidget* listWidget = ui->listWidget;
	listWidget->setViewMode(QListView::IconMode);
	listWidget->setIconSize(QSize(80, 80));
	listWidget->addItem(new QListWidgetItem(QIcon("/home/jast/FreeCAD/src/Mod/PartDesign/Gui/Resources/icons/bandsaw.jpg"), tr("Band Saw")));
	listWidget->addItem(new QListWidgetItem(QIcon("/home/jast/FreeCAD/src/Mod/PartDesign/Gui/Resources/icons/jigsaw.jpg"), tr("Jig Saw")));
	listWidget->addItem(new QListWidgetItem(QIcon("/home/jast/FreeCAD/src/Mod/PartDesign/Gui/Resources/icons/chospaw.jpg"), tr("Chop Saw")));
	listWidget->addItem(new QListWidgetItem(QIcon("/home/jast/FreeCAD/src/Mod/PartDesign/Gui/Resources/icons/drill.jpg"), tr("Drill")));
}

void ToolsConfiguration::InitializeTable()
{
	//ui->tableWidget->verticalHeader()->setVisible(false);
	ui->tableWidget->horizontalHeader()->setVisible(false);
}

void ToolsConfiguration::onChangedToolSelection(QListWidgetItem* w)
{
	currentListItem = w;
	int nRow = 0;
	ui->tableWidget->setColumnCount(2);
	ui->tableWidget->setRowCount(20);
	QTableWidgetItem *tableItem = new QTableWidgetItem("Type");
	InsertTableItem(tableItem, nRow, 0);
	tableItem = new QTableWidgetItem(w->text());
	InsertTableItem(tableItem, nRow++, 1);

	// Length 
	tableItem = new QTableWidgetItem("Length");
	InsertTableItem(tableItem, nRow, 0);
	tableItem = new QTableWidgetItem("1.0");
	InsertTableItem(tableItem, nRow++, 1);

	// Length 
	tableItem = new QTableWidgetItem("Length");
	InsertTableItem(tableItem, nRow, 0);
	tableItem = new QTableWidgetItem("1.0");
	InsertTableItem(tableItem, nRow++, 1);

	// Length 
	tableItem = new QTableWidgetItem("Length");
	InsertTableItem(tableItem, nRow, 0);
	tableItem = new QTableWidgetItem("1.0");
	InsertTableItem(tableItem, nRow++, 1);

	// Length 
	tableItem = new QTableWidgetItem("Length");
	InsertTableItem(tableItem, nRow, 0);
	tableItem = new QTableWidgetItem("1.0");
	InsertTableItem(tableItem, nRow++, 1);

	// Length 
	tableItem = new QTableWidgetItem("Length");
	InsertTableItem(tableItem, nRow, 0);
	tableItem = new QTableWidgetItem("1.0");
	InsertTableItem(tableItem, nRow++, 1);

	// Length 
	tableItem = new QTableWidgetItem("Length");
	InsertTableItem(tableItem, nRow, 0);
	tableItem = new QTableWidgetItem("1.0");
	InsertTableItem(tableItem, nRow++, 1);
}

void ToolsConfiguration::InsertTableItem(QTableWidgetItem* item, int rowIndex, int columnIndex)
{
	ui->tableWidget->setItem(rowIndex, columnIndex, item);
}

void ToolsConfiguration::onActionAddTool()
{
	if (currentListItem != nullptr)
	{
		QTreeWidgetItem* curTreeItem = new QTreeWidgetItem(QStringList(currentListItem->text()));
		ui->treeWidget->addTopLevelItem(curTreeItem);
	}
}

void ToolsConfiguration::setTabIndex(int idx)
{
	tabIdx = idx;
}

int ToolsConfiguration::getTabIndex() const
{
	return tabIdx;
}

#include "moc_ToolsConfiguration.cpp"
