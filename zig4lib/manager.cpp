#include "manager.h"

File::~File()
{
  // manager.remove(this);
}

void File::slotFileChanged()
{
  //qDebug() << "File::slotFileChanged()... ";
  emit fileChanged();
}

void intrusive_ptr_add_ref(File *p)
{
	++p->refCount;
}

void intrusive_ptr_release(File *p)
{
	if(--p->refCount == 0) {
		delete p;
	}
}

FileManager::FileManager(/* const FileChangedCallback& fccb */) /*: fccb(fccb) */
{
  // qDebug() << "FileManager::FileManager()";
  // connect(&watcher, SIGNAL(fileChanged(const QString&)), this, SLOT(slotFileChanged(const QString&)));
}

// void FileManager::watch(File *file)
// {
//   watcher.addPath(file->path); 
// }

// void FileManager::remove(File *file)
// {
//   // qDebug() << "FileManager::remove(): ";
//   if(map.get<FileTag>().erase(file) > 0) {
//     emit fileRemoved(file);
//   }
// }

// void FileManager::slotFileChanged(const QString& path) const
// {
//   //qDebug() << "FileManager::slotFileChanged()... ";
//   FileMap::iterator i(map.find(path));
//   if(i != map.end()) {
//     //qDebug() << "FileManager::slotFileChanged: " << path;
//     (*i)->slotFileChanged();
//   }
//   // fccb.doReact();
// }
