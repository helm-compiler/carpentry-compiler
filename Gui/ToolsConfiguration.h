#pragma once
# include <QDialog>
# include "ui_ToolsConfiguration.h"

class ToolsConfiguration :
	public QDialog
{
	Q_OBJECT

public:
	ToolsConfiguration(QWidget *parent = 0);
	~ToolsConfiguration();
	void InitializeTools();
	void InitializeTable();
	void setTabIndex(int idx);
	int getTabIndex() const;

public slots:
	void onChangedToolSelection(QListWidgetItem* w);
	void onActionAddTool();

private:
	void InsertTableItem(QTableWidgetItem* item, int rowIndex, int columnIndex);

public:
	Ui::ToolsConfig *ui;

private:
	QListWidgetItem *currentListItem;
	int tabIdx;
};

