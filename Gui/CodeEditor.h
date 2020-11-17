#ifndef CODEEDITOR_H
#define CODEEDITOR_H

#include <QPlainTextEdit>
#include <QObject>
#include "Highlighter.h"

class QPaintEvent;
class QResizeEvent;
class QSize;
class QWidget;

class LineNumberArea;


class CodeEditor : public QPlainTextEdit
{
	Q_OBJECT

public:
	CodeEditor(QWidget *parent = 0);
    ~CodeEditor() {};

	void lineNumberAreaPaintEvent(QPaintEvent *event);
	int lineNumberAreaWidth();
	void setTabIndex(int t);
	int getTabIndex() const;

protected:
	void resizeEvent(QResizeEvent *event) override;
	void highlightCurrentLine();

private Q_SLOTS:
	void updateLineNumberAreaWidth(int newBlockCount);
	void updateLineNumberArea(const QRect &, int);

private:
	QWidget *lineNumberArea;
	Highlighter *highlighter;
	int tabIndex;
};


class LineNumberArea : public QWidget
{
public:
	LineNumberArea(CodeEditor *editor) : QWidget(editor) {
		codeEditor = editor;
	}

	QSize sizeHint() const override {
		return QSize(codeEditor->lineNumberAreaWidth(), 0);
	}

protected:
	void paintEvent(QPaintEvent *event) override {
		codeEditor->lineNumberAreaPaintEvent(event);
	}

private:
	CodeEditor *codeEditor;
};


#endif